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
library(R2jags)
source("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV9.R")
source("/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstatV4.R")

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
# Analysis using MCMC

#Standardize the covariate
Spiders2$ColonySize.c <- Mystd(Spiders2$ColonySize)

# For the binomial coding in glmer, we need to specify 
# as response variable: cbind(Succes, Failure)
Spiders2$Failure <- Spiders2$Respondents - Spiders2$Attackers



#As explained in the description of this exercise, we want to
#fit the Bayesian version of the following model:
M2 <- glmer(cbind(Attackers, Failure) ~ Species + PreySizeClass + 
                                        ColonySize.c +
                                        PreySizeClass : ColonySize.c +
                            (1 | Colony) ,
           data = Spiders2,
           family = binomial)
############################################################


# Fit the binomial GLMM in JAGS
#1. Set up the JAGS.data object
X <- model.matrix(~ 1 + Species + PreySizeClass + ColonySize.c +
                    PreySizeClass : ColonySize.c, 
                  data = Spiders2)
K <- ncol(X) #Number of regression parameters

#Random effects colony
re  <- as.numeric(as.factor(Spiders2$Colony))
Numre <- length(unique(Spiders2$Colony))

#And put it all together
JAGS.data  <- list(Y      = Spiders2$Attackers,  #Response  
                   Trials = Spiders2$Respondents,#Trials  
                   X      = X,               #Intercept + covariates
                   K      = K,               #Number of betas       
                   N      = nrow(Spiders2),  #Sample size     
                   re     = re,              #Random effect   
                   Numre  = Numre)           #Number of random effects



#Step 2: Formulate JAGS modelling code
sink("BinomialGLMM.txt")
cat("
model{
    #1A. Priors regression parameters
    for (i in 1:K) { beta[i]  ~ dnorm(0, 0.0001) }  
              
    #1B. Priors for random intercepts
    for (i in 1:Numre) {a[i] ~ dnorm(0, tau) } 

    #1D. Priors for variance random intercepts
    sigma ~ dunif(0, 5)
    tau  <- 1 / (sigma * sigma)
    
    #2. Likelihood
    for (i in 1:N) {
       Y[i]          ~ dbin(Pi[i], Trials[i])  
       logit(Pi[i]) <- inprod(beta[], X[i,]) + a[re[i]] 
	}
                                     
    #3. Expectes values Y, variance Y, Pearson residuals
    for (i in 1:N) {    
       ExpY[i] <- Pi[i] * Trials[i] 
       VarY[i] <- Pi[i] * Trials[i] * (1 - Pi[i])
       E[i]    <- (Y[i]  - ExpY[i]) / sqrt(VarY[i])    
    }
       
    #Simulate data with mean/variance taken from the fitted model
    for (i in 1:N) {
       YNew[i] ~  dbin(Pi[i], Trials[i])                    
       ENew[i] <- (YNew[i] - ExpY[i]) / sqrt(VarY[i])       
       D[i]    <- pow(E[i], 2)                      
       DNew[i] <- pow(ENew[i], 2)   
     }               
    
    #Sum of squared Pearson residuals for observed and new data 
    Fit    <- sum(D[1:N])
    FitNew <- sum(DNew[1:N])
}
",fill = TRUE)
sink()
#####################################
#


#Set initial values
inits  <- function () {
  list(beta  = rnorm(K, 0, 0.01),     #betas
       a     = rnorm(Numre, 0, 0.1),  #random effects
       sigma = runif(1, 0, 5)         #sigma city
       )  }

#Step 4: Parameters to estimate
params <- c("beta", 
            "sigma", 
            "Fit", "FitNew", "E",
            "ExpY"
            )

#Step 5: Run JAGS
Bin1   <- jags(data        = JAGS.data,
                inits      = inits,
                parameters = params,
                model      = "BinomialGLMM.txt",
                n.thin     = 10, 
                n.chains   = 3,
                n.burnin   = 4000,
                n.iter     = 5000)

#And do some updating
Bin2 <- update(Bin1, n.iter = 10000, n.thin = 10)  
outB  <- Bin2$BUGSoutput

### Rhat and neff: Rhat measures variance across chains. Measure of 
### 1.0 indicates good mixing, but look at pictures. neff measures 
### autocorrelation within chains. max is n.iterations. increase iterations
### and or increase thinning to assure enough useful iterations


#Step 6: Present output 
print(outB, digits = 3)  


#Step 7: Assess mixing and check overdispersion
MyNames <- c(colnames(X), "sigma")
#Adjust this variable if extra parameters are added!

MyBUGSChains(outB, 
             c(uNames("beta", K), "sigma"),
              PanelNames = MyNames)
#We should do more iterations!

############################################
#Check for overdispersion / underdispersion
#Bayesian approach:
mean(outB$sims.list$Fit >  outB$sims.list$FitNew) 
#That is not ok...we have overdispersion!

#Plot the sum of squares for each iteration. 
Min <- min(outB$sims.list$Fit, outB$sims.list$FitNew)
Max <- max(outB$sims.list$Fit, outB$sims.list$FitNew)

xyplot(outB$sims.list$FitNew ~ outB$sims.list$Fit, 
       aspect = "iso", 
       col = 1, 
       pch = 16, 
       cex = 0.7,
       xlim = c(Min, Max), 
       ylim = c(Min, Max)) 
#Yes...overdispersion

#What does the frequentist way of thinking tell us?
Ep <- outB$mean$E
N <- nrow(Spiders2)
p <- K  + 1
sum(Ep^2) / (N - p)
# Yes...there is overdispersion...just a little bit!

# As in the frequentist approach, we now
# need to apply a detailed model validation
# and figure out why we have overdispersion
# You can copy and paste all the R code from the
# frequentist analysis.
E2 <- outB$mean$E
F2 <- outB$mean$ExpY


#Plot residuals vs fitted values
par(mfrow = c(1,1), cex.lab = 1.5, mar = c(5,5,2,2))
plot(x = F2,
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)
#No major problems.

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



#So..shall we just present the model as it is, 
#and mention in the discussion of a paper that 
#there is a tiny amount of overdispersion?
############################################################





############################################################
# Step 8: Numerical output
OUT1 <- MyBUGSOutput(outB, 
                     c(uNames("beta", K), "sigma"),
                     VarNames = MyNames)
print(OUT1, digits = 5)




#######################################################
# If there is higher overdispersion, then use either a 
# beta-binomial distribution
###############################################


##################################################
#So...what does the model tell us?
OUT1 <- MyBUGSOutput(outB, 
                     c(uNames("beta", K), "sigma"),
                     VarNames = MyNames)
print(OUT1, digits = 5)
#TASK: Write down the fitted model


# Visualize the model
# Create data on a grid, and the matching X matrix
MyData <- ddply(Spiders2, 
                .(Species, PreySizeClass), 
                summarize,
                ColonySize.c = seq(min(ColonySize.c), 
                                 max(ColonySize.c),
                                 length = 25))
head(MyData)

X <- model.matrix(~Species + PreySizeClass + ColonySize.c +
                                        PreySizeClass : ColonySize.c, 
                  data = MyData)

# Extract parameters
Beta.mcmc <- outB$sims.list$beta 

# Calculate for each MCMC iteration the predicted values on this grid.
eta.mcmc <- X %*% t(Beta.mcmc)
dim(X)
dim(Beta.mcmc)
dim(eta.mcmc)

Pi.mcmc <- exp(eta.mcmc) / (1 + exp(eta.mcmc))



# Why don't we take the 2.5% and 97.5% values
# at each of the artificial covariate values?
# And plot these instead of the 3,000 lines? 
# We wrote  a small support function that does this:
	
L <- GetCIs(Pi.mcmc)
MyData2 <- cbind(MyData,L)
head(MyData2)

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
p <- ggplot()
p <- p + geom_point(data = Spiders2, 
                    aes(y = Ratio.AR, x = ColonySize.c),
                    shape = 1, 
                    size = 1)
p <- p + xlab("ColonySize.c") + ylab("Probability of attacking")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData2, 
                   aes(x = ColonySize.c, y = mean), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData2, 
                     aes(x = ColonySize.c, 
                         ymax = up, 
                         ymin = lo ),
                     alpha = 0.2)
p <- p + facet_grid(PreySizeClass ~ Species, scales = "fixed")
p



#Or like this:
p <- ggplot()
p <- p + geom_point(data = Spiders2, 
                    aes(y = Ratio.AR, 
                        x = ColonySize.c,
                        colour = PreySizeClass),
                    shape = 1, 
                    size = 1)
p <- p + xlab("ColonySize.c") + ylab("Probability of attacking")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData2, 
                   aes(x = ColonySize.c, 
                       y = mean,
                       group = PreySizeClass,
                       colour = PreySizeClass))

p <- p + geom_ribbon(data = MyData2, 
                     aes(x = ColonySize.c, 
                         ymax = up, 
                         ymin = lo,
                         group = PreySizeClass,
                         fill = PreySizeClass),
                     alpha = 0.2)
p <- p + facet_grid(. ~ Species, scales = "fixed")
p


#TASK: Change these graph so that along the x-axis the original
#      scale is plotted

#Task: At which colony size value does the curve for 'large prey'
#      differ from the other two prey size classes?
#      How would you get a measure of uncertainty?


