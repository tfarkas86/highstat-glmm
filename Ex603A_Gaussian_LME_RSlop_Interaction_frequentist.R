#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#############################################################################
#Lilies exercise
#Set the working directory and import the data
setwd("/Users/Highstat/MaggieWindows/applicat/HighlandStatistics/webHighstatNew2_2010/CourseBaFrGLMM5/Data/")
LS <- read.table("Lilies.txt", 
                 header=TRUE)
#############################################################################
#Data from
#Beavers and lilies: selective herbivory and adaptive foraging
#behaviour
#Freshwater Biology (2014) 59, 224-232 doi:10.1111

#Based on pads collected from lakes where beavers are
#absent or where there is no evidence of grazing on
#N. alba, a strong positive relationship between
#the petiole diameter and midline distance of N. alba
#pads is to be expected!

names(LS)
#[1] "Site"              "Mid_line_Diameter" 
#[4] "Petiole_Diameter"  "Beaver"


str(LS)
#'data.frame':   356 obs. of  5 variables:
# $ Site             : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Mid_line_Diameter: num  118.5 87.6 95.9 98.9 99.1 ...
# $ Petiole_Diameter : num  6.47 5.74 5.76 5.48 6.53 4.3 5.93  ...
# $ Beaver           : Factor w/ 2 levels "N","Y": 1 1 1 1 1 1 1 1 1 1 ...

######################################################################
#Aim:
#Model Mid_line_Diameter  as a function of:
#  - Petiole_Diameter, Beaver absence/presence, interaction
#  - Site effect (random intercept + slope)

###################################################################
#Load packages and library files
library(lattice)
library(nlme)  
library(ggplot2)
source("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV9.R")
###################################################################



###################################################################
#Housekeeping
#Create a new variable with shorter names
LS$MD <- LS$Mid_line_Diameter
LS$PD <- LS$Petiole_Diameter

#Convert Site into a factor
LS$fSite <- factor(LS$Site)
###################################################################


###################################################################
str(LS)
table(LS$fSite)
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
# 15 20 19 34 21 20 20 24 31 12 20 20 20 20 20 20 20
###################################################################




######################################################################
#Data exploration

# OUTLIERS?
par(mfrow = c(2, 1), mar = c(3, 4, 1, 1))
dotchart(LS$MD, main= "Mid_line_Diameter")
dotchart(LS$PD, main= "Petiole_Diameter")

#No!
######################################################################
# RELATIONSHIPS?

plot(x = LS$PD, 
     y = LS$MD,
     pch = 16, cex = 0.5, 
     ylab = "Mid line Diameter (mm)",
     xlab = "Petiole Diameter (mm)")
#Perfect linear relationship!

#Site effect
boxplot(LS$MD ~ fSite,
        varwidth = TRUE, 
        data = LS,
        xlab     = "Site",
        ylab     = "Mid_line_Diameter")
abline(h = mean(LS$MD, na.rm = TRUE), 
       lty = 2)


#Is it possible that the MD - PD relationship
#changes by beaver absence/presence?
p <- ggplot()
p <- p + geom_point(data = LS, 
                    aes(y = MD, x = PD),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Petiole Diameter (mm)") + 
         ylab("Mid line Diameter (mm)")
p <- p + theme(text = element_text(size=15))
p <- p + geom_smooth(data = LS, 
                     aes(x = PD, y = MD), method = "lm")
p <- p + facet_grid( .~ Beaver, scales = "fixed")
p
#Hmmm....we have a sample size issue here!
#Maybe a simulation study to investigate the sample
#size effect may be wise!

#Is it possible that the MD - PD relationship
#changes by site?
p <- ggplot()
p <- p + geom_point(data = LS, 
                    aes(y = MD, x = PD),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Petiole Diameter (mm)") + 
         ylab("Mid line Diameter (mm)")
p <- p + theme(text = element_text(size=15))
p <- p + geom_smooth(data = LS, 
                     aes(x = PD, y = MD), method = "lm")
p <- p + facet_wrap( ~ fSite, ncol = 4, scales = "fixed")
p

#This one is nice too:
p <- ggplot()
p <- p + geom_point(data = LS, 
                    aes(y = MD, x = PD),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Petiole Diameter (mm)") + 
         ylab("Mid line Diameter (mm)")
p <- p + theme(text = element_text(size=15))
p <- p + geom_smooth(data = LS, 
                     aes(x = PD, 
                         y = MD,
                         group = fSite), 
                     method = "lm")
p <- p + facet_grid(. ~ Beaver, scales = "fixed")
p
#Hmmm...do we really want to fit this model?


#################################################################
#END OF DATA EXPLORATION
#################################################################

#Continue with frequentist analysis



################################################################
#Protocol for mixed modelling (3 steps).

##################################################################
#Step 1 of the protocol

#We fit the model:
# MD_ij    ~ N(mu_ij, sigma^2)
# E(MD_ij) = mu_ij
# mu_ij    = beta_1 + beta_2 * PD_ij + beta_3 * Beaver_ij +
#            beta4  * PD_ij * Beaver_ij + 
#            a_i + b_i * PD_ij

# a_i ~ N(0, sigma1_site^2)
# b_i ~ N(0, sigma2_site^2)

#That is a lot of sigmas!

#The implementation of such a model is given by:
M1 <- lme(MD ~ PD * Beaver,
          random =~ 1 + PD | fSite,  
          data = LS, 
          method = "REML")
summary(M1)
#Write out the model:

# MD_ij    ~ N(mu_ij, 12.17^2)
# E(MD_ij) = mu_ij

#Beavers not present:
# mu_ij    = -12.26 + 22.25 * PD_ij + a_i + b_i * PD_ij

#Beavers present
# mu_ij    = -12.26 + 7.11 + (22.25 - 3.15) * PD_ij + a_i + b_i * PD_ij

# But...the interaction is not significant!!!!



#Everything ok?
#Take residuals.....
#1.1 - Plot residuals versus fitted values
#1.2 - Plot residuals versus each covariate in the model
#1.3 - Plot residuals versus each covariate NOT in the model
#1.4 - Plot residuals versus time (if relevant)
#1.5 - Plot residuals versus spatial coordinates

#1.1
plot(M1, pch = 16, col = 1) ##OK

#Or
#Normalized residuals in lme are of the form:
# E / sqrt(variance of the residuals)

E1 <- resid(M1, type = "n")  #THIS IF FOR lme..NOT lme4
F1 <- fitted(M1)

par(mfrow = c(1, 1))
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals", 
     xlim = c(0, 220))
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1)

# Ok, no negative fitted values

#1.2
boxplot(E1 ~ LS$fSite, data = LS)
abline(0,0)

plot(x = LS$PD, y = E1)
abline(h = 0, lty = 2)


## Perfect, we do not see non-linear patterns

#1.3: Not relevant here

#1.4: Not relevant here

#1.5  Not relevant here

#Results model validation:
#Looks ok...no structure in residuals...go to step 2


##########################################################
#Step 2
#Is everything significant? 
#Approach C (likelihood ratio test)
M2A <- lme(MD ~ PD * Beaver,
          random =~ 1 + PD | fSite,  #This is just syntax!
          data = LS, 
          method = "ML")

M2B <- lme(MD ~ PD + Beaver,
          random =~ 1 + PD | fSite,
          data = LS, 
          method = "ML")
anova(M2A, M2B)

#The interaction is not significant. Drop it.



M3A <- lme(MD ~ PD + Beaver,
          random =~ 1 + PD | fSite,
          data = LS, 
          method = "ML")

M3B <- lme(MD ~ PD,
          random =~ 1 + PD | fSite,
          data = LS, 
          method = "ML")

M3C <- lme(MD ~ Beaver,
          random =~ 1 + PD | fSite,
          data = LS, 
          method = "ML")

anova(M3A, M3B)
#Beaver borderline significant.

#PD highly significant
anova(M3A, M3C)



#Let's have a look at the estimated parameters and 95% CI for
#Beaver presence.


M4 <- lme(MD ~ PD + Beaver,
          random =~ 1 + PD | fSite,
          data = LS, 
          method = "REML")
summary(M4)

#Task: Write down the estimated model

#Beavers not present:
#mu_i = -11.70 + 22.15 * PD_i + a_i + b_i * PD_i

#Beavers present:
#mu_i = -11.70 - 13.28 + 22.15 * PD_i + a_i + b_i * PD_i





########################################################################
#Step 3: Presents results (REML gives better estimates of the sigmas...
#        so use REML for final presentation)


#Just double check all the assumptions again (in case you 
#removed covariates)
plot(M4)
E4 <- resid(M4, type = "n")
#Repeat all model validation graphs from step 2
#.......








#######################################################################
#Sketch fitted values
#A. Specify covariate values for predictions
#B. Create X matrix with expand.grid
#C. Calculate predicted values
#D. Calculate standard errors (SE) for predicted values
#E. Plot predicted values
#F. Plot predicted values +/- 	1.96 * SE


#A:
range(LS$PD)
#[1] 1.57 9.93
MyData <- expand.grid(PD = seq(1.57, 9.93, length = 25),
                      Beaver = levels(LS$Beaver))

#Here is some fancy code to make a better MyData
#It doesn't extrapolate
library(plyr)
MyData <- ddply(LS, .(Beaver), summarize,
                PD = seq(min(PD), max(PD),length = 25))

#B. Create X matrix with expand.grid
X <- model.matrix(~ PD + Beaver, data = MyData)


#C. Calculate predicted values
#NewData$Pred <- predict(M4, NewData, level = 0)
#The level = 0 ensure that we fit the fixed effects
#Or:
MyData$Pred <- X %*% fixef(M4)



#D. Calculate standard errors (SE) for predicted values
#SE of fitted values are given by the square root of
#the diagonal elements of: X * cov(betas) * t(X)

MyData$SE <- sqrt(  diag(X %*% vcov(M4) %*% t(X))  )

#And using the Pred and SE values, we can calculate
#a 95% confidence interval
MyData$SeUp <- MyData$Pred + 1.96 * MyData$SE
MyData$SeLo <- MyData$Pred - 1.96 * MyData$SE


#E. Plot predicted values
MyData

p <- ggplot()
p <- p + geom_point(data = LS, 
                    aes(y = MD, x = PD),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Petiole Diameter (mm)") + 
         ylab("Mid line Diameter (mm)")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = PD, 
                       y = Pred), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = PD, 
                         ymax = SeUp, 
                         ymin = SeLo ),
                     alpha = 0.6)

p <- p + facet_grid(. ~ Beaver, scales = "fixed")
p

