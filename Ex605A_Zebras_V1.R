#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#############################################################################
#Zebra exercise
#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/NewDataSets/lme_ZebrasStress/")
Zebra <- read.csv("ZebraStress.csv", header=TRUE)
#############################################################################


# Data from:
# Diet quality in a wild grazer declines under the 
# threat of an ambush predator
# Barnier et al. (2015)
# Proc. R. Soc. B 281: 20140446.# http://dx.doi.org/10.1098/rspb.2014.0446


#From the abstract:
# Predators influence prey populations not only through predation 
# itself, but also indirectly through prompting changes in prey 
# behaviour. The behavioural adjustments of prey to predation 
# risk may carry nutritional costs, but this has seldom been 
# studied in the wild in large mammals. Here, we studied the# effects of an ambush predator, the African lion (Panthera leo), 
# on the diet quality of plains zebras (Equus quagga) in Hwange 
# National Park, Zimbabwe. We combined information on movements 
# of both prey and predators, using GPS data, and measurements 
# of faecal crude protein, an index of diet quality in the prey.# Zebras which had been in close proximity to lions had a lower 
# quality diet, showing that adjustments in behaviour when lions 
# are within short distance carry nutritional costs. The ultimate 
# fitness cost will depend on the frequency of predator–prey 
# encounters and on whether bottom-up or top-down forces are more 
# important in the prey population. Our finding is the first attempt# to our knowledge to assess nutritionally mediated risk effects 
# in a large mammalian prey species under the threat of an ambush 
# predator, and brings support to the hypothesis that the 
# behavioural effects of predation induce important risk effects 
# on prey populations.

#In short:
# 1. Zebras are tagged. Lions are tagged. So we know the distances
#    between them.
# 2. Zebra faecel samples are collected. Crude protein (CP) content 
#    is measured (CP in the diet is important for growth and 
#    reporduction in horses).
# 3. Setup of the experiment:
#       65 Facel samples from 30 individuals in 8 harems.
#       Multiple days.


# Response variable: CP
# Covariates: Sex (male, female). 
#             Distance to lions. Categorical (close, far) 
#                                This is minimum distance. 
#                                It is a predation index/risk

# Random effects:
# From the paper: 
# 1. We included a random effect on the intercept for 
#    data collected on the same day in the same group, 
#    thus accounting for the fact that samples coming 
#    from members of the same harem on the same date 
#    were correlated. 

#    Highstat: Coding required for the day-group identifier:
Zebra$fHaremDate  <- factor(paste(Zebra$Group, 
                                  Zebra$Date, 
                                  sep = "_"))
#and then (1 | fHaremDate) in lme4

# 2. We also included a random effect on the intercept 
#    for individuals nested in harems, to take account 
#    of repeated measures on individuals and groups.

#    Highstat: We think they used: (1 | fIndividual)
#    This would make sense to: (1 | fGroup / fIndividual)
#    but the small sample size may not allow it.

#    3. A random effect varying among individuals nested 
#       in harems was included on the coefficient relating 
#       distance to lions to CP in the faeces, to allow 
#       for different reactions to risk of predation between 
#       individuals in their groups,  

#     Highstat: We guess they used: (1 + fDistance | fGoup / fIndividual)
#     So..a random slope.

# So..the random error structre used in the paper is of the form:
# (1 | fHaremDate) + (1 + fDistance | fGroup / fIndividual) 
# We hope they have a large data set!


names(Zebra)
# [1] "Crude.Protein"     "Distance.to.lions" "Sex"   "Group"            
# [5] "Individual"        "Date"   


str(Zebra)

######################################################################
#Aim:
# Model CP as a function of Sex and distance to lions
# Take into account correlation.


###################################################################
#Load packages and library files

library(lattice)
library(lme4) 
library(ggplot2) 
source("/Users/Highstat/MaggieWindows/applicat/HighlandStatistics/webHighstatNew2_2010/CourseMCMCGLMM/RCode/HighstatLibV7.R")


###############################################
#Housekeeping: Define factors as factors
Zebra$fIndividual <- factor(Zebra$Individual)
Zebra$fGroup      <- factor(Zebra$Group)
Zebra$fHaremDate  <- factor(paste(Zebra$Group, 
                                  Zebra$Date, 
                                  sep = "_"))

Zebra$fDistance   <- factor(Zebra$Distance.to.lions,
                             levels = c("Far", "Close"))
Zebra$fSex        <- factor(Zebra$Sex)
Zebra$CP          <- Zebra$Crude.Protein



######################################################################
#Data exploration
#Balanced data?

table(Zebra$fIndividual) #Really?
table(Zebra$fHaremDate)  #Really?
table(Zebra$fDistance)   #Ok
table(Zebra$fSex)        #Ok


# Outliers?
MyVar <- c("CP")
Mydotplot(Zebra[, MyVar])
#No outliers!

# RELATIONSHIPS?
#Distance effect
boxplot(CP ~ fDistance,
        varwidth = TRUE, 
        data = Zebra,
        xlab     = "Distance to lions",
        ylab     = "CP")
#There seems tobe an effect!

boxplot(CP ~ fSex, data = Zebra) #Maybe

#And the random effects:
boxplot(CP ~ fHaremDate, data = Zebra)
boxplot(CP ~ fIndividual, data = Zebra)

#END OF DATA EXPLORATION
#################################################################


################################################################
#Continue with frequentist analysis


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


##################################################################
#Step 1 of the protocol

library(lme4)
# I think that the random slope stuff is too complicated for this
# data set.
M1 <- lmer(CP ~ fDistance + fSex + (1  | fIndividual) + (1 | fHaremDate),
           data = Zebra, 
           REML = TRUE)
summary(M1)

# Are distance and sex significant? 
# To answer this question, create 95% CIs,
# and notice that 0 is not on the 95% CI 
# interval for Distance, but the beta for 
# sex is. So distance is significant but sex 
# is not.

# Also note that the sigma for individual is 0.
# Should we drop it and simplify the random error 
# structure? Sounds sensible!


#Model validation.
#Everything ok?
#Take residuals.....
#1.1 - Plot residuals versus fitted values
#1.2 - Plot residuals versus each covariate in the model
#1.3 - Plot residuals versus each covariate NOT in the model
#1.4 - Plot residuals versus time (if relevant)
#1.5 - Plot residuals versus spatial coordinates

#1.1
E1 <- resid(M1)
F1 <- fitted(M1)

par(mfrow = c(1, 1))
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1)
#Some funny patterns for larger fitted values!
#No negative fitted values

#1.2
boxplot(E1 ~ fDistance, data = Zebra)
abline(h = 0, lty = 2)
#Good

boxplot(E1 ~ fSex, data = Zebra)
abline(h = 0, lty = 2)
#Ok-ish


#1.3: Not relevant here

#1.4: Not relevant here

#1.5  Not relevant here

#Results model validation:
#Looks ok...no structure in residuals...go to step 2


##########################################################
#Step 2
#Is everything significant? To answer
#this question we can use three different approaches in a
#frequentist approach:
#A. Use p-values of t-statistics (or z-statistics)
#B. Ue F test
#C. Use likelihood ratio test


#We have already done step A...and discovered that distance
#is significant, and sex is not.

#Approach C (likelihood ratio test)
#But the data set is small!!
M1a <- lmer(CP ~ fDistance + fSex + (1  | fIndividual) + (1 | fHaremDate),
           data = Zebra, 
           REML = FALSE)
drop1(M1a, test = "Chi")
#Same conclusion

########################################################################
#Step 3: Presents results (REML gives better estimates of the sigmas...
#        so use REML for final presentation).
#        So that is M1

#Understand the output:
summary(M1)

# Random effects:
 # Groups      Name        Variance  Std.Dev. 
 # fHaremDate  (Intercept) 4.908e-01 7.006e-01
 # fIndividual (Intercept) 8.363e-15 9.145e-08
 # Residual                2.091e-01 4.572e-01
# Number of obs: 65, groups:  fHaremDate, 34; fIndividual, 30

# Fixed effects:
               # Estimate Std. Error t value
# (Intercept)      5.6049     0.1655   33.87
# fDistanceClose  -0.7920     0.2993   -2.65
# fSexM           -0.2772     0.1629   -1.70

#For 'Far' and 'Females':
CP_ijk = 5.60

#For Far and Males:
CP_ijk = 5.60 - 0.277

#For Close and Females:
CP_ijk = 5.60 - 0.79 - 0.27

#The sigma_individual is nearly 0




###############################################################################
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
MyData <- expand.grid(fDistance = levels(Zebra$fDistance),
                      fSex      = levels(Zebra$fSex))

#B. Create X matrix with expand.grid
X <- model.matrix(~ fDistance + fSex, data = MyData)


#C. Calculate predicted values
MyData$Fit <- X %*% fixef(M1)



#D. Calculate standard errors (SE) for predicted values
#SE of fitted values are given by the square root of
#the diagonal elements of: X * cov(betas) * t(X)

MyData$SE <- sqrt(  diag(X %*%vcov(M1) %*% t(X))  )

#And using the Pred and SE values, we can calculate
#a 95% confidence interval
MyData$SeUp <- MyData$Fit + 1.96 * MyData$SE
MyData$SeLo <- MyData$Fit - 1.96 * MyData$SE


#E. Plot predicted values
MyData



p <- ggplot()
p <- p + xlab("Distance to lions") + ylab("CP")
p <- p + theme(text = element_text(size=15)) + theme_bw()

p <- p + geom_point(data = MyData, 
                    aes(x = fDistance, 
                        y = MyData$Fit, 
                        size = 6),    
                    col = ("red"))

p <- p + geom_errorbar(data = MyData,
                       aes(x = fDistance, 
                           ymax = SeUp, 
                           ymin = SeLo), 
                       width = 0.2,
                       col = "red")

p <- p + geom_point(data = Zebra, 
                    aes(x = fDistance, y = CP),
                    position = position_jitter(width = .02),
                    color = grey(0.3),
                    size = 2)
                           
p <- p + facet_grid(. ~ fSex, 
                    scales = "fixed")
p <- p + theme(legend.position="none") 
p <- p + theme(strip.text.y = element_text(size = 15, 
                                           colour = "black", 
                                           angle = 20),
               strip.text.x = element_text(size = 15, 
                                           colour = "black", 
                                           angle = 0)                            
                                           )
p
#

#What do you think of these analyses?

# It would be very nice to do a power analysis on these 
# data. If we were to repeat the experiment again, how 
# many observations would we need in order to find a 
# significant sex effect?




