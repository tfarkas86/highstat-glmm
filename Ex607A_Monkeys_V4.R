#    Highland Statistics Ltd.
#    www.highstat.com
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

#This exercise is based on:
# Evidence for varying social strategies across the day in 
# chacma baboons.
# Claudia Sick, Alecia J. Carter, Harry H. Marshall, Leslie A. Knapp, 
# Torben Dabelsteen and Guy Cowlishaw
# Biol. Lett. 2014 10, 20140249, published 9 July 2014
# Data where taken from the online supplemental material:
# http://rsbl.royalsocietypublishing.org/content/suppl/2014/07/08/rsbl.2014.0249.DC1.ht ml


#Description of experiment/data:
#1. Data were collected from two troops of chacma baboons
#2. Each troop was followed daily from dawn to dusk (ca 06.00–18.00 h). 
#3. 1 h long focal follows on all identifiable individuals 
#4. Recording all grooming and dominance interactions
#5. Unit of analysis:  grooming session
#6. Multiple observations from the same focal session
#7. Multiple observations on the same animal.
   

#Question 1: Subordinates should prefer to groom more dominant 
#            animals earlier in the day
#Question 2: Subordinates’ relative contributions to grooming sessions 
#            should be greater earlier in the day, especially with 
#            high-ranking partners


#What will you learn in this exercise?
# 1A. Reckognising where the correlation is (2 way nested and crossed data)
# 1B. Some simple ggplot2 coding
# 1C. lmer coding
# 1D. R^2 for a mixed model

# 2A. Binomial GLMM
# 3.  If this is part of an MCMC course, then a beta-GLMM will be coded



###########################################################################
#Set the working directory
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data")	
Monkeys <- read.table(file = "MonkeysV2.txt", 
                      header = TRUE,
                      dec = ".")
names(Monkeys)
str(Monkeys)
#####################################################################



#####################################################################
#Question 1: Subordinates should prefer to groom more dominant animals 
#            earlier in the day

#For this question we need to analyse only the Subordinate Grooms data:

Monkeys.sub <- Monkeys[Monkeys$SubordinateGrooms == "yes",]
MS          <- Monkeys.sub  #Is shorter


# The Sick et al. (2014) paper gives a description of the 
# experimental design, but at times it is difficult to 
# understand how exactly sampling took place.
# Here are the essential parts:


#1. Two monkey groups (large vs small). Coded as GroupSize
#2. There are multiple focal hours per day (coded as "FocalHour")
#3. In a focal hour we look for a groomer (coded as FocalGroomer)
#4. We also write down the id of the receiver monkey  (coded as Receiver)
#5. Record the time since sunrise, and the rank difference between 
#   focal groomer and receiver

#Aim: RankDifference = function(Time, Relatedness, GroupSize)

#Sensible interactions:
    #Time and relatedness  (used in paper)
    #Based on our common sense: Relatedness and groupsize  ?
#############################################################





###########################################################
# Understanding the dependency in the data
# Where is the correlation?
# Answer: See the pdf with the skecth of the design




# Sick et al use: 
# (1 | GroupSize /  FocalGroomer / FocalHour) + (1 | Receiver)
# The first one doesn't make sense (only 2 levels in GroupSize!).


# So...we better use:
# (1 | FocalGroomer / FocalHour) + (1 | Receiver)
# (1 | FocalGroomer) + ( 1 | FocalHour) + (1 | Receiver)


# FocalGroomer is groomer monkey identity. 
# FocalHour is the hour in which sampling took place



# First we need to figure out how things have been coded
table(MS$FocalGroomer)
table(MS$FocalHour)
length(table(MS$FocalGroomer))
length(table(MS$FocalHour))

# 59 monkeys divided over two groups
# 573 hours of sampling. Note: each sampling hour is
# uniquely coded! Each monkey is uniquely coded. 
# That makes life easy!


############################################
# Load lme4
library(lme4)



###########################################################################
#Response variable: 
# RankDifference  (difference between groomer and receiver)

#Covariates: Time           (time since sunrise)
#            Relatedness    (index specifying relatedness)
#            GroupSize      (small vs large)  
#
#            FocalGroomer   (monkey id who is doing the grooming)
#            FocalHour      (focal hour)
#            Receiver       (grooming partner identity)
#
#####################################################################




#####################################################################
#Housekeeping
library(lattice)
library(ggplot2)
source(file = "HighstatLibV9.R")
#####################################################################





######################################################################
#Data exploration

#Size of the data
dim(MS)

#Outliers
MyVar <- c("RankDifference", "Time", "Relatedness")
Mydotplot(Monkeys.sub[, MyVar])
# No outliers
# Rank difference is scaled between 0 and 1!
# What distribution should we use?



#Relationships
# Create a graph with Time along the x-axis and Rank differences
# along the y-axis
p <- ggplot(data = MS, aes(x = Time , y = RankDifference))
p <- p + geom_point()  #add points
p <- p + geom_smooth()  #Add a smoother
p <- p + facet_wrap(~ GroupSize)
p          # Non-linear pattern for the smaller group?


#Make a boxplot of RankDifference conditional on focal groomer 
p <- ggplot(MS, aes(y = RankDifference, x = FocalGroomer))
p <- p + geom_boxplot() 
p 
#Some groomers have consistently higher rank differences



#Collinearity
#Is there collinearity between Time and partner ID?
p <- ggplot(MS, aes(y = Time, x = Receiver))
p <- p + geom_boxplot() 
p    #Not really
######################################################################


######################################################################
# Frequentist analysis

#Based on the paper:
M1 <- lmer(RankDifference ~ Time * Relatedness + GroupSize +
                            (1 | FocalGroomer / FocalHour) + 
                            (1 | Receiver),
           data = MS)
summary(M1)


#Based on our common sense:
M2 <- lmer(RankDifference ~ Time * Relatedness + Relatedness * GroupSize +
                            (1 | FocalGroomer / FocalHour) + 
                            (1 | Receiver),
           data = Monkeys.sub)
summary(M2)


# RD_ijk = N(mu_ijk, sigma^2)

# E(RD_ijk) = mu_ijk
# mu_ijk = Intercept + Fixes Stuff + a_i + b_ij + c_k + eps_ijk

# i is Monkey index
# j is hour index 
# a_i ~  N(0, sigma_Monkey^2)
# b_ij ~ N(0, sigma_Hour^2)
# c_k  ~ N(0, sigma_Receiver^2)


#Due to the coding of the random effects in Excel, this is the same as:
M2B <- lmer(RankDifference ~ Time * Relatedness + Relatedness * GroupSize +
                            (1 | FocalGroomer) + 
                            (1 | FocalHour) + 
                            (1 | Receiver),
           data = Monkeys.sub)
summary(M2B)




##Model selection on the interactions
M2B.ml <- lmer(RankDifference ~ Time * Relatedness + 
                                Relatedness * GroupSize +
                            (1 | FocalGroomer) + 
                            (1 | FocalHour) + 
                            (1 | Receiver),
           data = Monkeys.sub,
           REML = FALSE)
drop1(M2B.ml, test = "Chi")
#Everything significant!
################################################################



################################################################
#Model validation
#Focus on model M2B

#Homogeneity
E2 <- resid(M2B)
F2 <- fitted(M2B)
plot(x = F2, 
     y = E2,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2)
#No clear heterogeneity..well..a little bit


#Independence
plot(y = E2,
     x = MS$Time,
     xlab = "Time",
     ylab = "Residuals")
abline(h = 0, lty = 2)
##Perfect..no non-linear pattern!
#Checking:
G1 <- gam(E2 ~ s(Time), data = MS)
summary(G1)
#No obvious time effect in the residuals


plot(y = E2,
     x = MS$Relatedness,
     xlab = "Relatedness",
     ylab = "Residuals")
abline(h = 0, lty = 2)
#Ok
#####################################################



#####################################################
#Model interpretation
# RankDifference_ijk ~ N(mu_ijk, sigma^2)
# E(RankDifference_ijk) = mu_ijk

# GroupSize = Large:
# mu_ijk = 0.63 - 0.008 * Time - 0.154 * Relatedness + 0.01 * Time * Relatedness

# GroupSize = Small:
# mu_ijk = 0.63 - 0.151 - 0.008 * Time + (0.225 - 0.154) * Relatedness + 
#          0.01 * Time * Relatedness


                            # Estimate Std. Error t value
# (Intercept)                 0.637661   0.037942  16.806
# Time                       -0.008282   0.002934  -2.823
# Relatedness                -0.154105   0.053522  -2.879
# GroupSizesmall             -0.151393   0.048497  -3.122
# Time:Relatedness            0.012807   0.006161   2.079
# Relatedness:GroupSizesmall  0.255740   0.044349   5.767





range(MS$Time)
range(MS$Relatedness)
MyData <- expand.grid(Time        = seq(0.266, 12.22, length = 15),
                      Relatedness = seq(0, 0.81, length = 15),
                      GroupSize   = levels(MS$GroupSize))

X          <- model.matrix(~Time * Relatedness + Relatedness * GroupSize, 
                           data = MyData)
Betas      <- fixef(M2)
MyData$Fit <- X %*% Betas
MyData$SE  <- sqrt(diag(X %*% vcov(M2) %*% t(X)))
MyData


#Skecth the fitted values:
for (i in 1: 300){
 p <- wireframe(Fit ~ Time + Relatedness, 
                data = MyData,
                group = GroupSize,
                #zlim = c(0,1),
                shade = TRUE,
                scales = list(arrows = FALSE),
                drape = TRUE, 
                colorkey = FALSE,
                screen = list(z = i, x = -60 - i /10))
 print(p)
}

 #Pick on value of i
 i <- 70
 p <- wireframe(Fit ~ Time + Relatedness, 
                data = MyData,
                group = GroupSize,
                zlim = c(0,1),
                shade = TRUE,
                scales = list(arrows = FALSE),
                drape = TRUE, 
                colorkey = FALSE,
                screen = list(z = i, x = -60 - i /5))
 print(p)
#Task: What does this mean in terms of ecology?
#############################################################



##############################################################
#Can we calculate an R^2?
#Answer: yes. 


#Paper: A general and simple method for obtaining R2 from#       generalized linear mixed-effects models
#       Shinichi Nakagawa and Holger Schielzeth
#       Methods in Ecology and Evolution. 2012
#       doi: 10.1111/j.2041-210x.2012.00261.x

#They show that to get an R^2 for a Gaussian mixed model, you
#need to carry out the following steps:

#1. Fit:  Y_ij = Fixed stuff + a_i + b_i + eps_ij
#2. Calculate variance of 'Fixed Stuff'
#3. Get the sigmas for the random effects a_i and b_i 
#4. Calculate two R^2s

#Variance explained by fixed effects AND random effects:

#           Variance Fixed Stuff + sigma_a^2 + sigma_b^2
# ---------------------------------------------------------------
#     Variance Fixed Stuff + sigma_a^2 + sigma_b^2 + sigma_eps^2


#Variance explained by fixed effects only:

#                         Variance Fixed Stuff 
# ----------------------------------------------------------------
#     Variance Fixed Stuff + sigma_a^2 + sigma_b^2 + sigma_eps^2


#We have three random intercepts.....just add a sigma_c^2


#To calculate these two R^2s:

X             <- model.matrix(~Time * Relatedness + 
                               Relatedness * GroupSize, 
                               data = Monkeys.sub)
Betas         <- fixef(M2)
Fixed         <- X %*% Betas
VarFixedStuff <- var(Fixed)

#And get the variance of the random effects
Sigma2.a   <- VarCorr(M2)$FocalGroomer[1]
Sigma2.b   <- VarCorr(M2)$'FocalHour:FocalGroomer'[1]
Sigma2.c   <- VarCorr(M2)$Receiver[1]
Sigma2.eps <- attr(VarCorr(M2), "sc")^2

#R2 based on covariates and random effects (conditional R2)
(VarFixedStuff + Sigma2.a + Sigma2.b + Sigma2.c) / 
(VarFixedStuff + Sigma2.a + Sigma2.b + Sigma2.b + Sigma2.eps)
#R2 = 81%

#R2 based on covariates only (marginal R2)
(VarFixedStuff) / (VarFixedStuff + Sigma2.a + Sigma2.b + Sigma2.c + Sigma2.eps)
#R2 = 5%  That is not much!!!
#So..the covariates are not doing much. 
#The random effects explain most of the 81%.



#Let us try to visualise why this is happening
#Plot fitted versus observed values
plot(x = fitted(M2), #Covariates + random effects
     y = MS$RankDifference,
     xlim = c(0,1),
     ylim = c(0,1))
#Good fit

#vs
plot(x = Fixed,      #Covariates only
     y = MS$RankDifference,
     xlim = c(0,1),
     ylim = c(0,1))


###############################################################################
#End of question 1
###############################################################################


