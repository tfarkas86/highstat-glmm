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
# Proc. R. Soc. B 281: 20140446.
# http://dx.doi.org/10.1098/rspb.2014.0446


#From the abstract:
# Predators influence prey populations not only through predation 
# itself, but also indirectly through prompting changes in prey 
# behaviour. The behavioural adjustments of prey to predation 
# risk may carry nutritional costs, but this has seldom been 
# studied in the wild in large mammals. Here, we studied the
# effects of an ambush predator, the African lion (Panthera leo), 
# on the diet quality of plains zebras (Equus quagga) in Hwange 
# National Park, Zimbabwe. We combined information on movements 
# of both prey and predators, using GPS data, and measurements 
# of faecal crude protein, an index of diet quality in the prey.
# Zebras which had been in close proximity to lions had a lower 
# quality diet, showing that adjustments in behaviour when lions 
# are within short distance carry nutritional costs. The ultimate 
# fitness cost will depend on the frequency of predator–prey 
# encounters and on whether bottom-up or top-down forces are more 
# important in the prey population. Our finding is the first attempt
# to our knowledge to assess nutritionally mediated risk effects 
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


#Continue with Bayesian analysis



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

#Reference model
M1 <- lmer(CP ~ fDistance + fSex + (1  | fIndividual) + (1 | fHaremDate),
           data = Zebra, 
           REML = TRUE)


#Start Bayesian code 
#A. Create the data for JAGS
X <- model.matrix(~ fDistance + fSex, data = Zebra)
K <- ncol(X)  #Number of columns
head(X)

#Random effects:
#We have two of these!
#This is the random effect fIndifidual
Ind     <- as.numeric(as.factor(Zebra$fIndividual))
Nre.Ind <- length(unique(Zebra$fIndividual))

#This is the random effect for the harem-day combination
Har     <- as.numeric(as.factor(Zebra$fHaremDate))
Nre.Har <- length(unique(Zebra$fHaremDate))

#So..now we have X * beta + (1 | fIndividual) + (1 | fHaremDate)


#Put it all together in a list
win.data <- list(Y       = Zebra$CP,    #Response
                 X       = X,           #Covariates
                 N       = nrow(Zebra), #Sample size
                 K       = K,           #Number of covariates
                 Ind     = Ind,         #random effect individual
                 Nre.Ind = Nre.Ind,     #number of random effects
                 Har     = Har,         #random effect harem-date
                 Nre.Har = Nre.Har)     #number of these random intercepts 
win.data

###################################################
# 2. JAGS modelling code

sink("zebra_lmm.txt")
cat("
model{
    #1A. Priors beta and sigma-epsilon
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}
    tau <- 1 / (sigma * sigma)
    sigma ~ dunif(0, 10)

    #1B. Priors random effects a1 (individual) and a2 (harem-date) 
    for (i in 1:Nre.Ind) { a1[i] ~ dnorm(0, tau_a1)}
    for (i in 1:Nre.Har) { a2[i] ~ dnorm(0, tau_a2)}
    
    #Matching sigmas for individual and harem-date
    #Use uniform priors
    tau_a1  <- 1 / (sigma_a1 * sigma_a1)
    tau_a2  <- 1 / (sigma_a2 * sigma_a2)
    sigma_a1 ~ dunif(0, 10)
    sigma_a2 ~ dunif(0, 10)

    #2. Likelihood
    for (i in 1:N) {
      Y[i]   ~ dnorm(mu[i], tau)   
      mu[i]  <- eta[i] + a1[Ind[i]] + a2[Har[i]]
      eta[i] <- inprod(beta[], X[i,]) 
 
      #3. Discrepancy measures 
      Exp[i] <- mu[i]      #Expected value of Y
      Var[i] <- sigma^2    #Variance of Y 
      
      #Residuals
      E[i]   <- (Y[i]  - Exp[i]) / sqrt(Var[i])    
      
     }          
}
",fill = TRUE)
sink()
#####################################
#

#3. Initial values & parameters to save
inits  <- function () {
  list(
    beta      = rnorm(K, 0, 0.1),        #betas
    sigma     = runif(1, 0, 10),         #sigma epsilon
    a1        = rnorm(Nre.Ind, 0, 0.1),  #random effect
    a2        = rnorm(Nre.Har, 0, 0.1),  #random effect
    sigma_a1  = runif(1, 0, 10),         #sigma individual
    sigma_a2  = runif(1, 0, 10))  }      #sigma harem-day


#What do we want to store?
params <- c("beta", "sigma", "E", "a1", "a2", 
            "sigma_a1", "sigma_a2", "mu", "eta")


#4. Start JAGS
J1 <- jags(data       = win.data,
           inits      = inits,
           parameters = params,
           model      = "zebra_lmm.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

J2 <- update(J1, n.iter = 10000, n.thin = 10)  
print(J2, digits = 3)

#5. Assess mixing
out <- J2$BUGSoutput


#Plot the chains and assess mixing
MyNames <- c(colnames(X), "Sigma", "Sigma individual", "Sigma group-day")

MyBUGSChains(out, 
             c(uNames("beta", K), "sigma", "sigma_a1", "sigma_a2"),
             PanelNames = MyNames)
 
#Auto-correlation in the chains?             
MyBUGSACF(out, 
          c(uNames("beta", K), "sigma", "sigma_a1", "sigma_a2"),
          PanelNames = MyNames)
          

#6. Present output
OUT1 <- MyBUGSOutput(out, 
                     c(uNames("beta", K), "sigma", "sigma_a1", "sigma_a2"),
                     VarNames = MyNames)
rownames(OUT1)[1:K] <- colnames(X)
print(OUT1, digits =5)

MyBUGSHist(out, 
           c(uNames("beta", K), "sigma", "sigma_a1", "sigma_a2"),
           PanelNames = MyNames)


OUT1
summary(M1)  #That is all good!
             #But note the difference in the sigma individual! 
######################################################




#Model validation
E   <- out$mean$E
Fit <- out$mean$mu 

plot(x = Fit, 
     y = E, 
     xlab = "Fitted values", 
     ylab = "Residuals")
abline(h = 0, lty = 2)


boxplot(E ~ Zebra$fSex)
abline(h = 0, lty = 2)

boxplot(E ~ Zebra$fDistance)
abline(h = 0, lty = 2)

##############################################################





#Step 2: Is everything significant?
rownames(OUT1)[1:K] <- colnames(X)
OUT1

MyBUGSHist(out, 
           c(uNames("beta", K), "sigma", "sigma_a1", "sigma_a2"),
           PanelNames = MyNames)

#No...not everything is significant!
#But keep the model as it is?
############################################################






########################################################
#Step 3: Explain the model
OUT1
#Write out the 4 equations....and sketch the fitted values




##########################################################
#Final task....sketch the fitted values with CI
# We have a large number of betas
# We can calculate the fitted values for each MCMC
# iteration and use these to calculate credible intervals:

betas <- out$sims.list$beta

#Create covariate values for prediction
MyData <- expand.grid(fDistance = levels(Zebra$fDistance),
                      fSex      = levels(Zebra$fSex))


#B. Create X matrix with expand.grid
Xp <- model.matrix(~ fDistance + fSex, data = MyData)


#Now calculate 6000 fitted values
mu.mcmc <- Xp %*% t(betas)

#Calculate the posterior mean and 95% credible intervals
# Why don't we take the 2.5% and 97.5% values
# at each of the artificial covariate values?
# And plot these instead of the 3,000 lines? 
# We wrote  a small support function that does this:
	
L <- GetCIs(mu.mcmc)
L
# GetCIs(): Support file to calculate the posterior
# mean, se and 2.5 and 97.5 quantiles for 
# fitted values. The input needs to be 
# an observation - by - MCMC iteration object.
 
# Each row in L is an artificial SDI value
# The first column is the posterior mean of 
# the MCMC iterations at that specific fitted 
# value.
# The second row is the SE, the third and fourth 
# are the 2.5 and 97.5% quartiles and are used for 
# the 95% credible interval. 

# Now plot the data and draw lines for the mean,
# and 2.5% and 97.5% values:

MyData2 <- cbind(MyData,L)
MyData2

#Each row in L is an artificial covariate 
#Now plot the data and draw lines for the mean,
#and 2.5% and 97.5% values:


library(ggplot2)
p <- ggplot()
p <- p + xlab("Distance to lions") + ylab("CP")
p <- p + theme(text = element_text(size=15)) + theme_bw()

p <- p + geom_point(data = MyData2, 
                    aes(x = fDistance, 
                        y = MyData2$mean, 
                        size = 6),    
                    col = ("red"))

p <- p + geom_errorbar(data = MyData2,
                       aes(x = fDistance, 
                           ymax = up, 
                           ymin = lo), 
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



