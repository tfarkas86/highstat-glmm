#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#############################################################################
#Lilies exercise
#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/")
LS <- read.table("Lilies.txt", 
                 header = TRUE)
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
#  - Petiole_Diameter
#  - Site effect

###################################################################
#Load packages and library files
library(lattice)
library(nlme)  
library(ggplot2)
source("HighstatLibV9.R")
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
p



#################################################################
#END OF DATA EXPLORATION
#################################################################

#Continue with frequentist analysis



################################################################
#Protocol for mixed modelling (3 steps).
#Based on Zuur et al. (2013).
#Beginner's Guide to GLM and GLMM with R

#1. Based on prior knowledge of the dependency
#   structure in the data, select the random
#   structure (i.e. the random effects) a priori.
#2. Fit the model and investigate which covariates
#   in the fixed part are significant or important.
#   Alternatively, omit model selection on the
#   covariates and adopt an information theoretic
#   approach (specify and compare 10-15 models).
#3. When the optimal model has been found, present
#   the numerical output and provide a graphic
#   representation of the model fit.

#In step 2 you may want to consider applying a limited amount of model
#selection...e.g. drop non-significant interactions.

##################################################################
#Step 1 of the protocol

#We first fit the random intercept model:
# MD_ij    ~ N(mu_ij, sigma^2)
# E(MD_ij) = mu_ij
# mu_ij    = beta_1 + beta_2 * PD_ij + a_i
# a_i ~ N(0, sigma_site^2)
 
M1 <- lme(MD ~ PD,
          random =~ 1 | fSite,  #This is just syntax!
          data = LS, 
          method = "REML")
summary(M1)

#But...what if the MD - PD relationship is changing per 
#site? Recall the ggplot2 figure!

# MD_ij    ~ N(mu_ij, sigma^2)
# E(MD_ij) = mu_ij
# mu_ij    = beta_1 + beta_2 * PD_ij + a_i + b_i * PD_ij
# a_i ~ N(0, sigma1_site^2)
# b_i ~ N(0, sigma2_site^2)

#That is a lot of sigmas!

#The implementation of such a model is given by:
M2 <- lme(MD ~ PD,
          random =~ 1 + PD | fSite,  #This is just syntax!
          data = LS, 
          method = "REML")
summary(M2)
#TASK: Write out the model

#Is the random intercept + slope model any better
#than the random intercept model?
AIC(M1, M2)
anova(M1, M2)

#But...what are we comparing here, and what is the null hypothesis?

# The random intercept model versus the 
# random intercept and slope model:
# MD_ij = alpha + fixed stuff + a_i + e_ij
# MD_ij = alpha + fixed stuff + a_i + b_i * PD_ij + e_ij

#a_i  ~ N(0, sigma1_site^2)
#b_i  ~ N(0, sigma2_site^2)
#e_ij ~ N(0, sigma^2)
#i is site index

anova(M1, M2)
#What are we testing here?
 #H0: sigma2_site = 0
 #H1: sigma2_site > 0

# This is called testing on the boundary as the sigma cannot be negative
# The p-value of the test is incorrect due to testing on the boundary.
# Approximate correction for testing on the boundary:
# 0.5 * (X2_1 + X2_2)

#run random intercept model           : call it M1
#run random intercept and slope model : call it M2
#run anova(M1, M2)
#Get test statistic.. L.ratio  and stick L.ratio in: 
0.5 * (  (1 - pchisq(8.880611, 1)) + (1 - pchisq(8.880611, 2))  )


#If you are going to use the random intercept and slope model....
#NOTE: INTRACLASS CORRELATION HAS A DIFFERENT DEFINITION
#This is why:
#y = X * beta + Z * b + eps
#Z = (ones  PD)
#cov(Y) = Z * cov(betas) * Z' + Sigma


#So..M2 is better...as judged by a p-value.
#Given the current climate that is not very favourable
#towards p-values, it would have been nicer to choose
#model M2 as starting model based on the underlying theory.
#Or based on model validation problems for model M1.



#Everything ok?
#Take residuals.....
#1.1 - Plot residuals versus fitted values
#1.2 - Plot residuals versus each covariate in the model
#1.3 - Plot residuals versus each covariate NOT in the model
#1.4 - Plot residuals versus time (if relevant)
#1.5 - Plot residuals versus spatial coordinates

#1.1
plot(M2, pch = 16, col = 1) ##OK

#Or
#Normalized residuals in lme are of the form:
# E / sqrt(variance of the residuals)

E2 <- resid(M2, type = "n")  #THIS IF FOR lme..NOT lme4
F2 <- fitted(M2)

par(mfrow = c(1, 1))
plot(x = F2, 
     y = E2,
     xlab = "Fitted values",
     ylab = "Residuals", 
     xlim = c(0, 220))
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1)

# Ok, no negative fitted values

#1.2
boxplot(E2 ~ LS$fSite, data = LS)
abline(0,0)

plot(x = LS$PD, y = E2)
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
M2A <- lme(MD ~ PD,
          random =~ 1 + PD | fSite,  #This is just syntax!
          data = LS, 
          method = "ML")

M2B <- lme(MD ~ 1,
          random =~ 1 + PD | fSite,
          data = LS, 
          method = "ML")
anova(M2A, M2B)


#Or with AIC (using ML!)
# AIC(M2A)
# AIC(M2B)
#etc....
#Finished step 2. 


########################################################################
#Step 3: Presents results (REML gives better estimates of the sigmas...
#        so use REML for final presentation)

M4 <- lme(MD ~ PD,
          random =~ 1 + PD | fSite,
          data = LS, 
          method = "REML")

#Just double check all the assumptions again (in case you 
#removed covariates)
plot(M4)
E4 <- resid(M4, type = "n")
#Repeat all model validation graphs from step 2
#.......


#Understanding the output:
# Mid_line_Diameter_ij = intercept + beta * Petiole_Diameter_ij  + 
#                        a_i + b_i * Petiole_Diameter_ij +  e_ij
# a_i ~ N(0, sigma1_Site^2)
# b_i ~ N(0, sigma2_Site^2)
# e_ij ~ N(0, sigma^2)
#   i = 1,....., 17 sites
#   j = 1, ..., n_i

summary(M4)
#We are also estimating a correlation between the
#random intercept and slope!

#Task: Write down the fitted model

# MD_i ~ N(mu_i, 12.22^2)
# E(MD_i)   = mu_i
# var(MD_i) = 12.22^2

# mu_i = -11.56 + 21.98 * PD_i + a_i + b_i * PD_i
       = -11.56 + a_i + ( 21.98 + b_i) * PD_i
# a_i ~ N(0, 5.83^2)
# b_i ~ N(0, 1.66^2)



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
MyData <- expand.grid(PD = seq(1.57, 9.93, length = 25))

#B. Create X matrix with expand.grid
X <- model.matrix(~ PD, data = MyData)


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

library(ggplot2)
p <- ggplot()
p <- p + geom_point(data = LS, 
                    aes(y = MD, x = PD),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Petiole Diameter (mm)") + 
         ylab("Mid line Diameter (mm)")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = PD, y = Pred), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = PD, 
                         ymax = SeUp, 
                         ymin = SeLo ),
                     alpha = 0.6)
p


##And if you want to have the 17 lines
LS$F0 <- fitted(M4, level = 0)  # alpha + beta * PD
LS$F1 <- fitted(M4, level = 1)  # alpha + beta * PD + 
                                # random intercept

p <- p + geom_line(data = LS, 
                   aes(x = PD, y = F1, group=fSite), 
                   colour = "black")
p
###########################################################








#########################################################
#B. Fit the same model using lmer from lme4 (need to be installed)
detach("package:nlme")  #nlme and lmer don't like each other
library(lme4)

#LM1 is equivalent of M1
LM1 <- lmer(MD ~ PD + (1 + PD | fSite),
            data = LS)

summary(LM1)

#How do you get p-values?
#Use a z-distribution
LM1Est.Param <- fixef(LM1)               #Get the betas
param.SD     <- sqrt(diag(vcov(LM1)))   #Get the SEs
pval         <- 2*pnorm(-abs(LM1Est.Param  / param.SD))  #Z distribution
cbind(LM1Est.Param, param.SD, pval)


#############################################################
#Start the protocol from scratch using lmer:

LM1 <- lmer(MD ~ PD + (1 + PD | fSite),
            data = LS)

#lmer graphical output: plot residuals vs fitted values
plot(fitted(LM1), resid(LM1))
abline(h = 0)


#Plot lme4 residuals
E1 <- resid(LM1)
boxplot(E1 ~ LS$fSite, data = LS)
abline(0,0)



#Looks ok...no structure in residuals...go to step 2
#Step 2: Is everything significant
LM2 <- lmer(MD ~ PD +
           (1 + PD | fSite), REML = FALSE,
          data = LS)

LM2A <- update(LM2, .~. -PD)
anova(LM2, LM2A)
#Or:
drop1(LM2, test = "Chi")
#Compare with lme results: Exactly the same

#Finished with model selection!
#Step 3
LM3 <- lmer(MD ~ PD + (1 + PD | fSite),
          REML = TRUE,
          data = LS)
summary(LM3)

#Model validation
plot(fitted(LM3), resid(LM3))
abline(h=0)

F3 <- fitted(LM3)
E3 <- resid(LM3)
plot(x = F3, y = E3)
plot(y = E3, x = LS$PD)
#End of protocol.





