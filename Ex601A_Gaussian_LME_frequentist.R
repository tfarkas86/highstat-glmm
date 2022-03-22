#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#################################################################
#Lilies exercise
#Set the working directory and import the data
#setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/")
LS <- read.table("Lilies.txt", 
                 header = TRUE,
                 dec = ".")
###############################################################
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
#Change the path of this:


library(lattice)
library(nlme)  #For lme
source("HighstatLibV9.R")

#Housekeeping
#Create a new variable with shorter names
LS$MD <- LS$Mid_line_Diameter
LS$PD <- LS$Petiole_Diameter

#Convert Site into a factor
LS$fSite <- factor(LS$Site)

str(LS)
table(LS$fSite)
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
# 15 20 19 34 21 20 20 24 31 12 20 20 20 20 20 20 20



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

# In step 2 you may want to consider applying a limited 
# amount of model selection...e.g. drop non-significant 
# interactions.

##################################################################
#Step 1 of the protocol
# Model:
#  MD_ij = beta1 + beta2 * PD_ij + a_i + eps_ij
#  a_i   ~ N(0, sigma_plot^2) 
#  eps_i ~ N(0, sigma^2) 
#
# Or: 
# MD_ij    ~ N(mu_ij, sigma^2)
# E(MD_ij) = mu_ij
# mu_ij    = Covariate stuff + a_i
#
# i = 1, .., 17
# j = 1, .., n_i

#
M1 <- lme(MD ~ PD,
          random =~ 1 | fSite,  #This is just syntax!
          data = LS, 
          method = "REML")

summary(M1)

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

E1 <- resid(M1, type = "n")  #THIS IS FOR lme..NOT lme4
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
abline(h = 0)

plot(x = LS$PD, y = E1)
abline(h = 0, lty = 2)

# Are there any patterns present?
# Do a GAM to check.
library(mgcv)
T1 <- gam(E1 ~ s(PD), data = LS)
summary(T1)
#No patterns present


## Perfect, we do not see non-linear patterns

#1.3: Not relevant here

#1.4: Not relevant here

#1.5  Not relevant here

#Results model validation:
#Looks ok...no structure in residuals...go to step 2


##########################################################
# Step 2
# Is everything significant? To answer
# this question we can use three different approaches in a
# frequentist approach:
# A. Use p-values of t-statistics (or z-statistics)
# B. Use F test
# C. Use likelihood ratio test (= drop1)

# Technicallity: In approaches A and B we MUST use REML estimation,
# in approach C we must use ML estimation

# If AIC...then use ML!
# ML/REML is only an issue for Gaussian mixed models. Not for GLMM

# Approach A:
M1 <- lme(MD ~ PD,
          random =~ 1 | fSite,
          data = LS, 
          method = "REML")

summary(M1) #And see what is significant based on p-values
##################################################
#End of approach A

#Approach B: F test

M2 <- lme(MD ~ PD,
          random =~ 1 | fSite,
          data = LS, 
          method = "REML")

anova(M2)  # But these p-values depend on the order of the 
           # covariates!
           # Ok in this case as there is only 1 covariate, and better perhaps 
           # if n:p is low, or you have orthogonal treatments
#End of approach B
##################################################

#Approach C (likelihood ratio test...= drop1)
M3 <- lme(MD ~ PD,
          random =~ 1 | fSite,
          data = LS, 
          method = "ML")


M3A <- update(M3, .~. -PD)
#M3A <- lme(MD ~ 1,
#           random =~ 1 | fSite,
#           data = LS, 
#           method = "ML")


# M3:  MD_ij = beta1 + beta2 * PD_ij + rest
# M3A: MD_ij = beta1 +                 rest
# H0: beta2 = 0


anova(M3, M3A)

#    Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#M3      1  4 2853.653 2869.153 -1422.827
#M3A     2  3 3647.505 3659.130 -1820.753 1 vs 2 795.8519  <.0001

#In paper: Likelihood ratio test: L = 795.85 (df = 1, p < 0.001)

# LRT not good if p:n is very high -- use Type 1 ANOVA instead


#Or with AIC (using ML!)
# AIC(M3)
# AIC(M3A)
#etc....
#Finished step 2. 


#############################################################
#Step 3: Presents results (REML gives better estimates of the 
#        sigmas...so use REML for final presentation)

M4 <- lme(MD ~ PD,
          random =~ 1 | fSite,
          data = LS, 
          method = "REML")

#Just double check all the assumptions again (in case you 
#removed covariates)
plot(M4)
E4 <- resid(M4, type = "n")
#Repeat all model validation graphs from step 2
#.......

boxplot(E4 ~ fSite, data=LS)
plot(E4 ~ PD, data=LS)


#Understanding the output:
# Mid_line_Diameter_ij = intercept + beta * Petiole_Diameter_ij  + 
#                        a_i + e_ij
# a_i  ~ N(0, sigma_Site^2)
# e_ij ~ N(0, sigma^2)
#   i = 1,....., 17 sites
#   j = 1, ..., n_i

summary(M4)

 #Linear mixed-effects model fit by REML
 #Data: LS
 #      AIC      BIC    logLik
 # 2850.317 2865.794 -1421.158

#Random effects:
# Formula: ~1 | fSiteNew
#        (Intercept) Residual
#StdDev:    8.512763 12.48395

#Fixed effects: MD ~ PD
#                Value Std.Error  DF  t-value p-value
#(Intercept) -13.40424 3.1032460 338 -4.31943       0
#PD           22.25932 0.3950169 338 56.35031       0
# Correlation:
#   (Intr)
#PD -0.714

#Standardized Within-Group Residuals:
#       Min         Q1        Med         Q3        Max
#-2.5925867 -0.6036509 -0.0073010  0.5763095  3.1548742

#Number of Observations: 356
#Number of Groups: 17

# a_i ~  N(0, 8.51^2)
# e_ij ~ N(0, 12.48^2)
#
# Intraclass correlation
#
#                  8.512763^2
#          -------------------------
#          8.512763^2 + 12.48395^2
#
# (8.512763 ^ 2)/( 8.512763^ 2 + 12.48395^2)     =  0.3173982




##############################################################
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
head(X)

#C. Calculate predicted values
#NewData$Pred <- predict(M4, NewData, level = 0)
#The level = 0 ensure that we fit the fixed effects
#Or:
MyData$Pred <- X %*% fixef(M4)  # = X * beta



#D. Calculate standard errors (SE) for predicted values
#   SE of fitted values are given by the square root of
#   the diagonal elements of: X * cov(betas) * t(X)  
#   Take this for granted!

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
p <- p + theme(text = element_text(size=15))
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
#A. What is what in nlme?
#   In this block of code we calculate the
#   fitted values and residuals ourselves.

# It turns out that nlme fitted values include both the fixed and random 
# effects. 

#   Extract the random effects:
RandomEffectSite <- ranef(M4)$'(Intercept)'
Site.Numerical   <- as.numeric(as.factor(LS$fSite))
RandomEffectSite[Site.Numerical]

#Y = X*beta + Z*b + eps
OUT <- cbind(LS$MD,
             fitted(M4,level = 0), #X * beta
             fitted(M4,level = 1), #X * beta + random effect
             RandomEffectSite[Site.Numerical],   #b
             resid(M4, level = 0),    #Z*b + eps
             resid(M4, level = 1),    #eps
             LS$fSite,                #Site number
             resid(M4, type = "n"),   #eps / sigma
             resid(M4)                #eps
             )

colnames(OUT) <- c("Y",
                 "X*beta",
                 "X*beta + Z*b",
                 "b",
                 "Z*b + eps",
                 "eps",
                 "Site",
                 "eps/sigma",
                 "eps")

head(OUT,10)  #First 10 rows of Z


#By defaul level = 1 in fitted() and resid()
#fitted(M4)[1:5]

# fitted(M4, level = 0)     : X * beta
# fitted(M4, level = 1)     : X * beta + Z * b    level = 1 is default
# ranef(M4)                 : b
# ranef(M4, level = 1)      : b
# resid(M4, level = 0)      : Z*b + eps
# resid(M4, level = 1)      : eps    level = 1 is default
# resid(M4)                 : eps
# resid(M4, type = "n")     : eps / sigma
# resid(M4, type = "n")     : resid(M7)/M7$sigma




#########################################################
#B. Fit the same model using lmer from lme4 (need to be installed)
detach("package:nlme")  #nlme and lmer don't like each other
library(lme4)

#LM1 is equivalent of M1
LM1 <- lmer(MD ~ PD + (1 | fSite),
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

LM1 <- lmer(MD ~ PD + (1 | fSite),
            data = LS)

#lmer graphical output: plot residuals vs fitted values
plot(fitted(LM1), resid(LM1))
abline(h = 0)

quartz()  #Mac command for a new window
#win.graph() on windows
plot(M1, col = 1)  #This is the lme graph


#Same question again....what is in resid and fitted????
#We have Y = X * beta + Z * b + epsilon
#Text below shows:
#   resid() is eps
#   fitted() is X * beta + Z * b


#In lme:
plot(M1)  #lme: These residuals are divided by sigma


#Plot lme4 residuals
E1 <- resid(LM1)
boxplot(E1 ~ LS$fSite, data = LS)
abline(0,0)



#Looks ok...no structure in residuals...go to step 2
#Step 2: Is everything significant
LM2 <- lmer(MD ~ PD+
           (1 | fSite), REML = FALSE,
          data = LS)

LM2A <- update(LM2, .~. -PD)
anova(LM2, LM2A)
#Or:
drop1(LM2, test = "Chi")
#Compare with lme results: Exactly the same

#Finished with model selection!

#Step 3
LM3 <- lmer(MD ~ PD + (1 | fSite),
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




#######################################################################
#The code below will not be discussed during the
#course...consider it as homework.

#What is what in lmer?
#Y = X*beta + Z*b + eps
LM3 <- lmer(MD ~ PD + (1 | fSite),
          REML = FALSE, data = LS)

#Get lmer fitted values
F.LM3 <- fitted(LM3)

#Calculate lmer fitted values manual
beta         <- fixef(LM3)
X            <- model.matrix(LM3)
F.LM3.manual <- X %*% beta

#Comparee lmer fitted values and manual calculates residuals
cbind(F.LM3, F.LM3.manual)  #These are different
#Why? The random effects are included in fitted.
#To see this, add the random effects to our F.LM3.manual:

#These are the 17 random effects
RE.LM3 <- ranef(LM3)$fSite$'(Intercept)'

#Blow them up to 356 rows
AllRE <- RE.LM3[as.numeric(LS$fSite)]

#And compare everything
Z.LM3 <- cbind(F.LM3, F.LM3.manual, AllRE,
               LS$Site, F.LM3.manual + AllRE )
colnames(Z.LM3) <- c("fitted","fit manual",
                     "Site","RE",
                     "fit manual + RE")

head(Z.LM3)
#Conclusions:
#fitted(lmer object) equal Xb + Zb


#Then the same holds for the residuals
resid(LM3)[1:5]  #Seems that these are the raw residuals (Y-Xb -Zb)
ResRawManual <- LS$MD - F.LM3.manual - AllRE

#Yes....resid(lmer object) for Gaussian model is Y - Xb - Zb
#Makes sense as there are no variance structures.



Z.LM3 <- cbind(F.LM3, F.LM3.manual, AllRE,
               LS$Site, F.LM3.manual + AllRE,
               resid(LM3), ResRawManual)
colnames(Z.LM3) <- c("fitted","fit manual",
                     "Site","RE",
                     "fit manual + RE",
                     "resid", "eps")

head(Z.LM3)


#If you want to divide things by sigma, use:
Sigmas <- c(sqrt(VarCorr(LM3)$fSite[1]), 
            attr(VarCorr(LM3), "sc")) 
Sigmas  



Z.LM3 <- cbind(F.LM3, F.LM3.manual, AllRE,
               LS$Site, F.LM3.manual + AllRE,
               resid(LM3), ResRawManual,
               resid(LM3)/Sigmas[2])
colnames(Z.LM3) <- c("fitted","fit manual",
                     "Site","RE",
                     "fit manual + RE",
                     "resid", "eps",
                     "resid/sigma")

head(Z.LM3)





