#    Highland Statistics Ltd.
#    www.highstat.com
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

# This exercise is based on:
# Differences in group size and the extent of individual 
# participation in group hunting may contribute to differential 
# prey-size use among social spiders
# Gyan Harwood and Leticia Aviles
# Biol. Lett. 2013 9, 20130621, published 27 November 2013


# From the abstract:
# We have previously shown that the range of prey sizes 
# captured by co-occurring species of group-hunting social 
# spiders correlates positively with their level of sociality. 
# Here, we show that this pattern is probably caused by 
# differences among species in colony size and the extent 
# to which individuals participate in group hunting.
# .....


# Description of experiment/data:
# 1. Three spider species examined (they co-occur in a forest in 
#    Sao Paulo, Brazil). 
# 2. Number of colonies per species:
#      16 for A. baeza 
#      29 for A. jabaquara
#      32 for A. dubiosus
# 3. Each colony is repeatedly presented with a prey (small, medium , large)
# 4. They measured:
#       - Total number of spiders in a colony
#       - Total number of spiders that emerged and surrounded the 
#         prey (= respondents)
#       - Number of attackers.
# 
#    The paper uses the following response variables:
#      A. 'respondents - attackers'. They called this the 
#         'attackers-to-observers differential
#      B. Total number of attackers from each hunt.
#    They used a Gaussian mixed model and a Poisson GLMM for these 
#    variables. Both are incorrect.

# Covariates in the model:
#   Colony size (log10)
#   Species 
#   Prey-size class
#   Plus all interactions
#   Plus a random effect colony identity

# Task: Model the probability of attackers as a function of 
#       the covariates. This is a binomial GLMM:

# Attackers_ij      ~ Bin(Pi_ij, Trials_ij)
# E(Attackers_ij)   = Pi_ij * Trials_ij
# var(Attackers_ij) = Pi_ij * Trials_ij * (1 - Pi_ij)
# Trials_ij is number o responders

# logit(Pi_ij) = Covariate stuff + a_i
# a_i ~ N(0, sigma^2_Colony)
# a_i is the random intercept for colony
######################################################




###########################################################################
#Set the working directory and load the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/")	
Spiders <- read.table(file = "SpidersGroupHuntingV2.txt", 
                      header = TRUE,
                      dec = ".")
names(Spiders)
str(Spiders)
##########################################################################



##########################################################################
# Load packages and support files
library(lattice)
library(ggplot2)
library(lme4)
library(plyr)
source("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV9.R")
##########################################################################



##########################################################################
#Housekeeping
Spiders$Attackers   <- Spiders$TotalNumberofSpidersThatAttackedPrey
Spiders$Respondents <- Spiders$TotalNumberRespondents
Spiders$Ratio.AR    <- Spiders$Attackers / Spiders$Respondents 
Spiders$ColonySize  <- Spiders$Population
Spiders$ColonySize.c <- MyStd(Spiders$ColonySize)


# So... Respondents is the number of trials
# It is the Trials_ij in the equation above.


# So we have Attackers times success out of Respondents 
# independent trials. If you don't want to use a binomial
# distribution, then model the ratio Attackers / Respondents
# with a beta distribution. See Zuur et al. (2012) for
# a detailed example.
############################################################



############################################################
#Data exploration
# In a binomial distirbution the number of trials needs to
# be larger than 0

Spiders$Respondents

# On a few occasions we have 0 respondents. These have to be removed from
# the analsyis!

Spiders2 <- subset(Spiders, Respondents > 0)
dim(Spiders)
dim(Spiders2)


#Size of the data
dim(Spiders2)

#How many observations per level of the factors
table(Spiders2$Species)  
table(Spiders2$Colony)   
table(Spiders2$PreySizeClass)
#Thats seems all ok
######################################################


#Outliers
MyVar <- c("Attackers", "Respondents", "Ratio.AR", "ColonySize")
Mydotplot(Spiders2[, MyVar])
# Quite often they all attack. 
# No outliers

#Collinearity
boxplot(ColonySize ~ PreySizeClass,
        data = Spiders2)
#No trouble        

boxplot(ColonySize ~ Species,
        data = Spiders2)
        
#No trouble        


##############
#Relationships
plot(y = Spiders2$Ratio.AR,
     x = Spiders2$ColonySize)

boxplot(Ratio.AR ~ Species,
        data = Spiders2)

boxplot(Ratio.AR ~ PreySizeClass,
        data = Spiders2)
#This doesn't make you happy!


#Interations
p <- ggplot(data = Spiders2, aes(x = ColonySize , y = Ratio.AR))
p <- p + geom_point()  
p <- p + facet_grid(PreySizeClass ~ Species)
p   
# We have enough data for some interactions!
# But do we have biological questions that guide us towards
# including interactions?
######################################################################





######################################################################
# Analysis

# For the binomial coding in glmer, we need to specify 
# as response variable: cbind(Succes, Failure)
Spiders2$Failure <- Spiders2$Respondents - Spiders2$Attackers

M1 <- glmer(cbind(Attackers, Failure) ~ Species + PreySizeClass + 
                                        ColonySize.c +
                                        Species : PreySizeClass + 
                                        Species : ColonySize.c +
                                        PreySizeClass : ColonySize.c +
                            (1 | Colony) ,
           data = Spiders2,
           family = binomial)
summary(M1)
# Hmmm...some of these SEs look suspiciously similar!


# Check for overdispersion
E1 <- resid(M1, test = "pearson")
N <- nrow(Spiders2)
p <- length(fixef(M1)) + 1
sum(E1^2) / (N - p)

# This is a bit on the high side!
# And some of the SEs look suspicious!


# Fit the same model in glmmADMB and compare results!
library(glmmADMB)
M1.admb <- glmmadmb(cbind(Attackers, Failure) ~ Species + PreySizeClass + 
                                        ColonySize.c +
                                        Species : PreySizeClass + 
                                        Species : ColonySize.c +
                                        PreySizeClass : ColonySize.c +
                                        (1 | Colony) ,
           data = Spiders2,
           family = "binomial")

summary(M1.admb)$coefficients #glmmADMB
summary(M1)$coefficients      #lme4

# Better print them next to each other
# Betas of admb and glmer:
cbind(summary(M1.admb)$coefficients[,1],
      summary(M1)$coefficients[,1])

# SEs of admb and glmer
cbind(summary(M1.admb)$coefficients[,2],
      summary(M1)$coefficients[,2])


# I think we can go on with glmer results
# Apply a likelihood ratio test:
drop1(M1, test = "Chi")




# If you do a backwards selection using the AIC you end up with:
M2 <- glmer(cbind(Attackers, Failure) ~ Species + PreySizeClass + 
                                        ColonySize.c +
                                        PreySizeClass : ColonySize.c +
                            (1 | Colony) ,
           data = Spiders2,
           family = binomial)
drop1(M2, test = "Chi")


E2 <- resid(M2, test = "pearson")
N <- nrow(Spiders2)
p <- length(fixef(M2)) + 1
sum(E2^2) / (N - p)

# There is still a little bit overdispersion. Let's do
# model validation.
######################################################




######################################################
# Model validation

# Get the fitted values
F2 <- fitted(M2)

# Plot residuals vs fitted values
par(mfrow = c(1,1), cex.lab = 1.5, mar = c(5,5,2,2))
plot(x = F2,
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)
#No major problems.

# Plot residuals vs each covariate
par(mfrow = c(1,1), cex.lab = 1.5, mar = c(5,5,2,2))
plot(x = Spiders2$ColonySize.c,
     y = E2,
     xlab = "Colony size",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)
#Ok

par(mfrow = c(1,1), cex.lab = 1.5, mar = c(5,5,2,2))
boxplot(E2 ~ Species,
     xlab = "Species",
     ylab = "Pearson residuals",
     data = Spiders2)
abline(h = 0, lty = 2)
#Perfect

par(mfrow = c(1,1), cex.lab = 1.5, mar = c(5,5,2,2))
boxplot(E2 ~ PreySizeClass,
     xlab = "PreySizeClass",
     ylab = "Pearson residuals",
     data = Spiders2)
abline(h = 0, lty = 2)
#Perfect

par(mfrow = c(1,1), cex.lab = 1.5, mar = c(5,5,2,2))
boxplot(E2 ~ Colony,
     xlab = "Colony",
     ylab = "Pearson residuals",
     data = Spiders2)
abline(h = 0, lty = 2)
# Hmmm..it would be nice to have the spatial locations
# of the colonies! Some of the boxplots next to each other
# in the graph are similar in behaviour. Are they also
# next to each other in the field?



# Why did they divide PreySize into classes? Why not using
# PreySizemm?
library(mgcv)
G1 <- gam(E2 ~ s(PreySizemm), data = Spiders2)
summary(G1)
plot(G1)  # Ok....there is not a nonlinear pattern 
          # present in the residuals

# So..shall we just present the model as it is, 
# and mention in the discussion of a paper that 
# there is a tiny amount of overdispersion?




#######################################################
# If there is higher overdispersion, then use either a 
# beta-binomial distribution, or the observation level 
# random intercept. The later one goes like this:
Spiders2$Eps <- 1:nrow(Spiders2)
M3 <- glmer(cbind(Attackers, Failure) ~ Species + PreySizeClass + 
                                        ColonySize.c +
                                        PreySizeClass : ColonySize.c +
                                        (1 | Colony) + (1 | Eps),
           data = Spiders2,
           family = binomial)
summary(M3)
drop1(M3, test = "Chi")
#We will not work this out any further.
###############################################




###############################################
#So...what does M2 tell us?
summary(M2)
#TASK: Write down the fitted model

# Attackers_ij ~ Bin(Pi_ij, N_ij)
# E(Attackers_ij) = Pi_ij * N_ij

# logit(Pi_ij) = eta_ij + a_i

# Baseline species, baseline prey size
# eta_ij = 1.82 + 0.60 * ColonySize.c





#And let's sketch the model fit.
M2 <- glmer(cbind(Attackers, Failure) ~ Species + PreySizeClass + 
                                        ColonySize.c +
                                        PreySizeClass : ColonySize.c +
                            (1 | Colony) ,
           data = Spiders2,
           family = binomial)



#Create data on grid, and the matching X matrix
range(Spiders2$ColonySize.c)
MyData <- expand.grid(Species       = levels(Spiders2$Species),
                      ColonySize.c  = seq(-0.87, 3.90, length = 25),
                      PreySizeClass = levels(Spiders2$PreySizeClass))

#Or better (to avoid extrapolation):
MyData <- ddply(Spiders, 
                .(Species, PreySizeClass), summarize,
                ColonySize.c = seq(min(ColonySize.c), 
                                 max(ColonySize.c),
                                 length = 25))
head(MyData)



X <- model.matrix(~Species + PreySizeClass + ColonySize.c +
                   PreySizeClass : ColonySize.c, 
                  data = MyData)
#Extract parameters and parameter covariance matrix
betas    <- fixef(M2)


#Calculate the fitted values in the predictor scale
MyData$eta <- X %*% betas
MyData$Pi  <- exp(MyData$eta) / (1 + exp(MyData$eta))

#Calculate the SEs on the scale of the predictor function
MyData$se    <- sqrt(diag(X %*% vcov(M2) %*% t(X)))
MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$se) / 
                (1 + exp(MyData$eta  + 1.96 *MyData$se))
MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$se) / 
                (1 + exp(MyData$eta  - 1.96 *MyData$se))

head(MyData)

#And make a picture
p <- ggplot()
p <- p + geom_point(data = Spiders2, 
                    aes(y = Ratio.AR, x = ColonySize.c),
                    shape = 1, 
                    size = 1)
p <- p + xlab("ColonySize.c") + ylab("Probability of attacking")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = ColonySize.c, y = Pi), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = ColonySize.c, 
                         ymax = SeUp, 
                         ymin = SeLo ),
                     alpha = 0.2)
p <- p + facet_grid(PreySizeClass ~ Species, scales = "fixed")
p



#Or like this:
p <- ggplot()
p <- p + geom_point(data = Spiders2, 
                    aes(y = Ratio.AR, x = ColonySize.c),
                    shape = 1, 
                    size = 1)
p <- p + xlab("ColonySize.c") + ylab("Probability of attacking")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = ColonySize.c, 
                       y = Pi,
                       group = PreySizeClass,
                       colour = PreySizeClass))

p <- p + geom_ribbon(data = MyData, 
                     aes(x = ColonySize.c, 
                         ymax = SeUp, 
                         ymin = SeLo,
                         group = PreySizeClass,
                         fill = PreySizeClass),
                     alpha = 0.2)
p <- p + facet_grid(. ~ Species, scales = "fixed")
p






#####################################
#Wanne have an R2?
source(file = "/Users/Highstat/applicat/HighlandStatistics/Courses/Data/R2glmms.R")
r.squared.merMod(M2)  
     # Class   Family  Link  Marginal Conditional      AIC
# 1 glmerMod binomial logit 0.1173288   0.2923799 670.9089

# So...we have an R2 of 12% (covariate effect only)
#################################################




