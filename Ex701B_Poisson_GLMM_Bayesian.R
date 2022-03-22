#    Introduction to Introduction to MCMC, 
#    Linear Mixed Effects models and GLMM with R
#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
#Rudy turnstone example
#Reference:
#Fuller RA, Bearhop S, Metcalfe NB, Piersma T (2013) 
#The effect of group size on vigilance in Ruddy Turnstones 
#Arenaria interpres varies with foraging habitat. 
#Ibis 155: 246-257

#4.1 Group size effect on vigilance in ruddy turnstones
#The data used in this chapter were presented in Fuller et al. (2013) 
#who investigated the effect of group size on vigilance in ruddy 
#turnstones Arenaria interpres in relation to the quality of 
#foraging habitat (in terms of food abundance and the threat of 
#their own predation) along a broadly linear 40-km stretch of 
#coastline in northeastern England.
#The ruddy turnstone is a short-legged shorebird species that 
#gets its name from its habit of turning over stones or clearing 
#away debris when it searches for intertidal food such as crustaceans, 
#worms, molluscs, and insects. 
#In their paper, Fuller et al. (2013) compared vigilance in habitat 
#types that differed greatly in prey abundance and proximity to 
#cover from which predators could launch surprise attacks.
#################################################################




#################################################################
#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Books/BGS/GAMM/Data/TurnStoneData")
TS <- read.table(file = "TurnstoneDataV2.txt", 
                 header = TRUE)

names(TS)
# [1] "BirdID"          "FlockID"         "HeadUps"         "FlockSize"       "TideState"      
# [6] "NumPecks"        "TotalAggression" "TimeOfDay"       "TimeHighTide"    "DaysStartWinter"
#[11] "Temperature"    
###################################################################




###################################################################
#Load packages and library files
library(lattice)  
library(lme4)
source(file = "/Users/Highstat/applicat/HighlandStatistics/Courses/FilesOnlineFollowUpRegressionGLMGAM/FinalExercises/HighstatLibV6.R")  
source("/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstat.R")
library(R2jags)
##################################################################


##################################################################
#Housekeeping
TS$fFlockID   <- factor(TS$FlockID)
TS$fTideState <- factor(TS$TideState)        
TS$TideStat01 <- as.numeric(TS$TideState)        
##################################################################



##################################################################
#Data exploration
#Outliers
MyVar <- c("HeadUps", "FlockSize", "NumPecks", "TimeOfDay",
           "TimeHighTide", "DaysStartWinter", "Temperature")
Mydotplot(TS[,MyVar])
#No extreme outliers

#Collinearity
MyVar <- c("FlockSize", "NumPecks", "TimeOfDay",
           "TimeHighTide", "DaysStartWinter", "Temperature",
           "TideStat01")

Mypairs(TS[,MyVar])
corvif(TS[,MyVar])

boxplot(FlockSize ~ TideState, data = TS)
boxplot(NumPecks ~ TideState, data = TS)
boxplot(TimeOfDay.Proportion. ~ TideState, data = TS)
boxplot(MinsSinceHighTide ~ TideState, data = TS)
boxplot(DaysSinceStartOfWinter ~ TideState, data = TS)
boxplot(Temperature ~ TideState, data = TS)

#Some trouble going on

#Relationships
MyVar <- c("FlockSize", "NumPecks", "TimeOfDay",
           "TimeHighTide", "DaysStartWinter", "Temperature")

Myxyplot(TS,MyVar,"HeadUps", MyYlab = "Head ups", MyXlab ="Covariates" )
#Some weak non-linear patterns! It is a GAMM chapter..:-)

#Balanced design?
par(mar = c(5,5,2,2))
plot(table(TS$fFlockID), 
     type = "h",
     xlab = "Flock identity",
     ylab = "Number of observations per flock",
     cex.lab = 1.5)


#Zero inflation
sum(TS$HeadUps == 0) / nrow(TS)



#####################################################
#Start MCMC analysis
#In this example we will show:
# 1. How to run a Poisson GLMM in JAGS


#Use standardized covariates
MyStd <- function(x) {(x - mean(x)) / sd(x)}
                     
TS$FlockSize.std    <- MyStd(TS$FlockSize)
TS$NumPecks.std     <- MyStd(TS$NumPecks)
TS$TimeOfDay.std    <- MyStd(TS$TimeOfDay)
TS$TimeHighTide.std <- MyStd(TS$TimeHighTide)
TS$Temperature.std  <- MyStd(TS$Temperature)
                     
#By design of the study: Use flock id as random effect
#The reference model is:

detach(package:gamm4)
detach(package:mgcv)
detach(package:nlme)
library(lme4)
M1 <- glmer(HeadUps ~ fTideState + 
                     FlockSize.std + 
                     NumPecks.std + 
                     TimeOfDay.std +
                     TimeHighTide.std  + 
                     Temperature.std + (1 | fFlockID),
           data = TS, family = poisson)

summary(M1)

##################################################
#Now fit the Poisson GLMM in JAGS

#Step 1: Fit the model.
#Use MCMC to estimate the parameters                     
#1. Bundle data
X <- model.matrix(~ fTideState + 
                     FlockSize.std + 
                     NumPecks.std + 
                     TimeOfDay.std +
                     TimeHighTide.std  + 
                     Temperature.std , data = TS)
   
                   
K <- ncol(X)  #Number of columns
head(X)

#Random effects:
Flock <- as.numeric(as.factor(TS$fFlockID))
Nre   <- length(unique(TS$fFlockID))


#The code below is copied an pasted from the owls
#Nest was changed into flock
win.data <- list(Y     = TS$HeadUps, 
                 X     = X,
                 N     = nrow(TS),
                 K     = K,
                 Flock = Flock,
                 Nre   = Nre)
win.data

###################################################
# 2. JAGS modelling code

#This code was copied from the coral reef example. We only changed Site into Nest.
sink("turnstonemmJAGS.txt")
cat("
model{
    #1A. Priors beta
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}

    #1B. Priors random effects and sigma_Flock
    # uniform or gamma sigmas (norm) and thetas (nb) in mixed models lead to poor mixing
    # half-Cauchy is meant to be better, but JAGS doesn't have half-Cauchy
    # if you divide two normally distributed random variables, you approximate
    # a half-Cauchy, so that's what is happening below for the assignment of 
    # sigma priors is to define two noramlly distributed random variables and
    # and divide them to make a half-Cauchy distributed sigma.

    for (i in 1:Nre) { a[i] ~ dnorm(0, tau_Flock)}
    
    num          ~ dnorm(0, 0.0016)      #<----half-Cauchy(25)
    denom        ~ dnorm(0, 1)           #<----half-Cauchy(25)
    sigma_Flock <- abs(num / denom)      #<----half-Cauchy(25)
    tau_Flock   <- 1 / (sigma_Flock * sigma_Flock)


    #2. Likelihood
    for (i in 1:N) {
      Y[i]        ~ dpois(mu[i])  
      log(mu[i]) <- eta[i]  
      eta[i]     <- inprod(beta[], X[i,]) + a[Flock[i]] 
  
      #3. New stuff: Discrepancy measures 
      Exp[i] <- mu[i] 
      Var[i] <- mu[i]
      E[i]   <- (Y[i]  - Exp[i]) / sqrt(Var[i])    
      #Pearson residuals
     
      YNew[i] ~  dpois(mu[i])                      
      #Simulated data with mean/variance taken from the 
      #fitted model
      
      ENew[i] <- (YNew[i] - Exp[i]) / sqrt(Var[i]) 
      #Normalized residual for predicted data
      
      D[i]    <- pow(E[i], 2)                      
      DNew[i] <- pow(ENew[i], 2)   
     }          
     Fit         <- sum(D[1:N])    #Sum of squared residuals  
     FitNew      <- sum(DNew[1:N]) #Sum of squared predicted 
                                   #residuals
}
",fill = TRUE)
sink()
#####################################
#

#A few more lines of code we copied from the coral
#reef example. Again..we replaced Site by Nest
#3. Initial values & parameters to save
inits  <- function () {
  list(
    beta  = rnorm(K, 0, 0.01),
    a     = rnorm(Nre, 0, 0.01),
    num   = rnorm(1, 0, 25),   #<-- half-Cauchy(25)
    denom = rnorm(1, 0, 1)  )  }


params <- c("beta", "E", "a", 
            "sigma_Flock", "mu", 
            "Fit", "FitNew")


#Don't forget to change the file name!
#4. Start JAGS
J1 <- jags(data       = win.data,
           inits      = inits,
           parameters = params,
           model      = "turnstonemmJAGS.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

J2 <- update(J1, n.iter = 10000, n.thin = 10)  
out <- J2$BUGSoutput

print(J2, digits = 3)  #Takes a few minutes
#Write down DIC at this point


#5. Assess mixing
MyBUGSChains(out, c(uNames("beta", K), "sigma_Flock"))
MyBUGSACF(out, c(uNames("beta", K), "sigma_Flock"))
#There is now some auto-correlation!


#6. Present output
OUT1 <- MyBUGSOutput(out, c(uNames("beta", K), "sigma_Flock"))
print(OUT1, digits = 5)

MyBUGSHist(out, c(uNames("beta", K), "sigma_Flock"))
#Some parameters are not significant!


OUT1
summary(M1)  #That is all good!
######################################################



######################################################
# Assess goodness of fit: Bayesion p-value:
# We have a large number of FitNew and Fit values. 
# Count how often one is bigger than the other:

mean(out$sims.list$Fit >  out$sims.list$FitNew)
#He....that is strange!

E   <- out$mean$E
sum(E^2) / (nrow(TS) - 8)
#Strange...so that overdispersion of 1.2 is doing
#some damage here!

#You can also plot the sum of squares for each iteration. 
Min <- min(out$sims.list$Fit, out$sims.list$FitNew)
Max <- max(out$sims.list$Fit, out$sims.list$FitNew)

xyplot(out$sims.list$FitNew ~ out$sims.list$Fit, 
       aspect = "iso", 
       col = 1, 
       pch = 16, 
       cex = 0.7,
       xlim = c(Min, Max), 
       ylim = c(Min, Max)) 


#################################################################



             
#Model validation
E   <- out$mean$E
Fit <- out$mean$mu #(includes the random effects)



#Plot residuals versus fitted values
plot(x = Fit, 
     y = E, 
     xlab = "Fitted values", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)


#Plot residuals vs each covariate
plot(x = TS$FlockSize, 
     y = E, 
     xlab = "FlockSize", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(x = TS$NumPecks, 
     y = E, 
     xlab = "NumPecks", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(x = TS$TimeOfDay, 
     y = E, 
     xlab = "TimeOfDay", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(x = TS$TimeHighTide, 
     y = E, 
     xlab = "TimeHighTide", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(x = TS$Temperature, 
     y = E, 
     xlab = "Temperature", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(x = TS$fTideState, 
     y = E, 
     xlab = "fTideState", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

#Step 3:
#Explain what it all means.

#Write down the fitted equation.



