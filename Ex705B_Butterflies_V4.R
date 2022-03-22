#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#######################################################################

#These data were taken from:
# A resource-based conservation approach for an 
# endangered ecotone species: the Ilex 
# Hairstreak (Satyrium ilicis) in 
# Flanders (north Belgium)

# Dirk Maes, Ilf Jacobs, Natascha Segers, Wouter Vanreusel 
# Toon Van Daele, Guy Laurijssens, Hans Van Dyck

########################################################

setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data")
BF <- read.table("ButterfliesNEggs_V3.txt", header = TRUE)
names(BF)
str(BF)


###########################################################################
#Load support files and packages
source(file = "HighstatLibV9.R")
source("/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstatV3.R")

library(lme4)
library(lattice)
library(R2jags)
library(rgl)
###########################################################################



# Response: Number of eggs.
# Covariates: We have a whole bunch of covariates, 
#             and their names are self-explanatory.
# We have multiple observation per 'landscape'. That 
# is like a site.

# The variables "Distance2Edge" and "NearestTallOak" are
# expressed in meters.
###########################################################################




###########################################################################
#Data exploration

#Outliers
MyVar <- c("NEggs", "LandscapeID", "NLowBranches",
          "TreeHeight", "HerbCover", "Shelter", 
          "BuckthornAbundance", "BrambleAbundance", "SmallOakAbundance",   
          "Distance2Edge", "NearestTallOak", "X",  "Y" )
Mydotplot(BF[,MyVar])
#Nothing too alien


#Zero inflation!
table(BF[,"NEggs"])


###################################
#Define factors as factors
BF$fShelter     <- factor(BF$Shelter)
BF$fLandscapeID <- factor(BF$LandscapeID)

#How many observations per cluster, and per shelter level?
table(BF$LandscapeID)
table(BF$fShelter)

#Is there a landscape effect?
boxplot(NEggs ~ LandscapeID, data = BF)


#Collinearity
MyVar <- c("NLowBranches",
          "TreeHeight", "HerbCover",
          "BuckthornAbundance", "BrambleAbundance", "SmallOakAbundance",   
          "Distance2Edge", "NearestTallOak", "X",  "Y" )

corvif(BF[,MyVar])
Mypairs(BF[,MyVar])


#Spatial position of the sampling locations
#X and Y are in Belgium coordinates
xyplot(Y ~ X,
       data = BF,
       aspect = "iso")
###################################################




###################################################
#MCMC analysis
MyStd <- function(x) { (x - mean(x)) / sd(x)}

BF$NLowBranchesc  <- MyStd(BF$NLowBranches)
BF$TreeHeightc    <- MyStd(BF$TreeHeight)
BF$HerbCoverc     <- MyStd(BF$HerbCover)
BF$BuckthornAbundancec <- MyStd(BF$BuckthornAbundance)
BF$BrambleAbundancec   <- MyStd(BF$BrambleAbundance)
BF$SmallOakAbundancec  <- MyStd(BF$SmallOakAbundance)
BF$Distance2Edgec      <- MyStd(BF$Distance2Edge)
BF$NearestTallOakc     <- MyStd(BF$NearestTallOak)
                    

# The aim of this exercise is to 
# get further familiar with JAGS code for Poisson GLMM and
# NB GLMM. We decided to select two covariates that were 
# the most important in the frequentist analysis. And we 
# will only use these covariates in the analysis. It makes 
# computing time shorter and it simplifies the model 
# interpretation part. For a full analysis you have to 
# include all covariates, or deal in another way with the 
# covariates (e.g. IT approach, forward selection, back ward 
# selection, etc.)


#########################################
#Fit a Poisson GLMM in JAGS

#1. Set up the win.data object
#Count part and binary parts of the model
X <- model.matrix(~ 1 + NLowBranchesc + Distance2Edgec, data = BF)
K <- ncol(X) #Number of regression parameters count part

#Random effects landscape
reLands  <- as.numeric(as.factor(BF$LandscapeID))
NumLands <- length(unique(BF$LandscapeID))

#And put it all together
win.data  <- list(Y        = BF$NEggs,  #Response    
                  X        = X,         #Intercept + covariates
                  K        = K,         #Number of betas       
                  N        = nrow(BF),  #Sample size     
                  reLands  = reLands,   #Random effect landscape
                  NumLands = NumLands)  #Random effect landscape 
win.data


#Step 2: Formulate JAGS modelling code
sink("PoissonGLMM.txt")
cat("
model{
    #1A. Priors regression parameters count part and binary part
    for (i in 1:K) { beta[i]  ~ dnorm(0, 0.0001) }  
              
    #1B. Priors for random intercepts landscape
    for (i in 1:NumLands) {a[i] ~ dnorm(0, tau.La) } 

    #1D. Priors for variances for random intercepts
    #    In case of non-mixing, change to half-Cauchy(25) or 
    #    half-Cauchy(16)
    sigma.La ~ dunif(0, 5)
    tau.La  <- 1 / (sigma.La * sigma.La)
    
    ###################################
    #2. Likelihood
    for (i in 1:N) {
       Y[i] ~ dpois(mu[i])
    	
       #log-link 
       log(mu[i]) <- eta[i]
       eta[i]     <- inprod(beta[], X[i,]) + a[reLands[i]]  

                                          
       #3. Discrepancy measures: 
       #   Expectes values, variance, Pearson residuals
       expY[i] <- mu[i] 
       varY[i] <- mu[i] 
       PRes[i] <- (Y[i] - expY[i]) / sqrt(varY[i])   
       
       #Simulate data from a ZIP  
       YNew[i]   ~   dpois(mu[i]) 
       PResNew[i] <- (YNew[i] - expY[i]) / sqrt(varY[i])
       D[i]       <- pow(PRes[i], 2)
       DNew[i]    <- pow(PResNew[i], 2)
    }     
     Fit         <- sum(D[1:N])
     FitNew      <- sum(DNew[1:N])
}
",fill = TRUE)
sink()
#####################################
#


#Inits function
#Set initial values
inits  <- function () {
  list(beta      = rnorm(K, 0, 0.01),        #betas
       a         = rnorm(NumLands, 0, 0.1),  #random effects
       sigma.La  = runif(1, 0, 5)            #sigma landscape
       )  }


#Step 4: Parameters to estimate
params <- c("beta", 
            "sigma.La", 
            "Fit", "FitNew", "PRes",
            "expY"
            )

#Step 5: Run JAGS
Pois1   <- jags(data        = win.data,
                inits      = inits,
                parameters = params,
                model      = "PoissonGLMM.txt",
                n.thin     = 10, 
                n.chains   = 3,
                n.burnin   = 4000,
                n.iter     = 5000)

#And do some updating
Pois2  <- update(Pois1, n.iter = 10000, n.thin = 10)  
#Pois3  <- update(Pois2, n.iter = 50000, n.thin = 10)  
#save(Pois3, file = "BF_Pois3.RData")

#It would be better to do at least another 100K iterations, 
#but this means that not all participants would be able 
#to finish the analysis.

#To reload results from an earlier session, use:
#load(file = "BF_Pois3.RData")
#recompile(ZIP3)
#load.module("dic")

out <- Pois2$BUGSoutput
#If you deciced to use Pois3, then adjust the code
#above.



#Step 6: Present output 
print(out, digits = 3)  


#Step 7: Assess mixing and check overdispersion
MyNames <- c(colnames(X), "sigma clutch")
#Adjust this variable if extra parameters are added!

MyBUGSChains(out, 
             c(uNames("beta", K), "sigma.La"),
              PanelNames = MyNames)
# We need more iterations here! But that requires half
# an hour computing time


############################################
#Check for overdispersion / underdispersion
#Bayesian approach:
mean(out$sims.list$Fit >  out$sims.list$FitNew) 
#That is not ok...we have overdispersion!

#Plot the sum of squares for each iteration. 
Min <- min(out$sims.list$Fit, out$sims.list$FitNew)
Max <- max(out$sims.list$Fit, out$sims.list$FitNew)

xyplot(out$sims.list$FitNew ~ out$sims.list$Fit, 
       aspect = "iso", 
       col = 1, 
       pch = 16, 
       cex = 0.7,
       xlim = c(Min, Max), 
       ylim = c(Min, Max)) 
#Yes...overdispersion

#What does the frequentist way of thinking tell us?
Ep <- out$mean$PRes
N <- nrow(BF)
p <- K  + 1
sum(Ep^2) / (N - p)



#As in the frequentist approach, we now
#need to apply a detailed model validation
#and figure out why we have overdispersion


#Step 8: Numerical output
OUT1 <- MyBUGSOutput(out, 
                     c(uNames("beta", K), "sigma.La"),
                     VarNames = MyNames)
print(OUT1, digits = 5)




#Your homework:
#1. Deal with model selection. We choose these covariates
#   based on some initial analyses.


#Task:
#Apply a model validation


#To understand why we have ovedispersion,
#it can also be handy to plot the model fit



#Calculate the fitted values.
#1. Get the betas from the MCMC iterations
#2. Define a grid of covariate values
#3. Calculate the predicted values on this 
#   grid for each MCMC iteration.
#4. Calculate the 95% credible intervals.
#5. Plot the whole thing


#1. Get the betas and gammas
beta.mcmc  <- out$sims.list$beta 
dim(beta.mcmc)

#Using these we can puzzle together the pieces.
#Note that we are not adding the random effects
#to the count part.

#2. Define a grid of covariate values
range(BF$NLowBranchesc)
range(BF$Distance2Edgec)

MyData <- expand.grid(NLowBranchesc  = seq(-2.15, 2.02, length = 25),
                      Distance2Edgec = seq(-0.96, 1.55, length = 25))
                       
#B. Create X matrix with expand.grid
X <- model.matrix(~ 1 + NLowBranchesc + Distance2Edgec, data = MyData)


#C. Calculate the predicted MCMC values
#   for mu
eta.mcmc  <- X %*% t(beta.mcmc) 
mu.mcmc   <- exp(eta)



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

ExpY <- MyData2[,"mean"]  #Expected values Y. For a ZIP these will change. 



#So..these are the two variables that we use
#to calculate our grid for:
BD.25   <- seq(-2.15, 2.02, length = 25)
Prc.25  <- seq(-0.96, 1.55, length = 25)
        
#And we convert the vector with expected values into
#a 25 by 25 matrix

ExpY.2d <- matrix(ExpY, nrow = length(Prc.25), ncol = length(BD.25))


#Now plot the jittered data again
library(rgl)
set.seed(12345)
BF$Covariate1.jitter <- jitter(BF$NLowBranchesc, amount = 0.1)
BF$Covariate2.jitter <- jitter(BF$Distance2Edgec, amount = 0.1)

plot3d(x = BF$Covariate1.jitter ,
       y = BF$Covariate2.jitter,
       z = BF$NEggs,
       type = "p",
       size = 10,
       lit = FALSE,
       xlab = "NLowBranchesc",
       ylab = "Distance2Edgec",
       zlab = "Number of eggs"
       )

#Add the surface for the fitted Poisson values
surface3d(BD.25, Prc.25, ExpY.2d, 
          alpha = 0.6, 
          front = "lines", 
          back = "lines", 
          color = "black"
          )

###################################################################




#Now fit a NB GLMM in JAGS

#Reference model
library(glmmADMB)
BF$fLandscapeID <- factor(BF$LandscapeID)
M2 <- glmmadmb(NEggs ~ NLowBranchesc +  Distance2Edgec +
                    (1 | fLandscapeID),
          family = "poisson",
          zeroInfl = FALSE,
          data = BF)
summary(M2)



#1. Set up the wind.data object
#Count part and binary parts of the model
X <- model.matrix(~ 1 + NLowBranchesc + Distance2Edgec, data = BF)
K <- ncol(X) #Number of regression parameters count part

#Random effects landscape
reLands  <- as.numeric(as.factor(BF$LandscapeID))
NumLands <- length(unique(BF$LandscapeID))

#And put it all together
win.data  <- list(Y        = BF$NEggs,  #Response    
                  X        = X,         #Intercept + covariates
                  K        = K,         #Number of betas       
                  N        = nrow(BF),  #Sample size     
                  reLands  = reLands,   #Random effect landscape
                  NumLands = NumLands)  #Random effect landscape 
win.data


#Step 2: Formulate JAGS modelling code
sink("NBGLMM.txt")
cat("
model{
    #1A. Priors regression parameters count part and binary part
    for (i in 1:K) { beta[i]  ~ dnorm(0, 0.0001) }  
              
    #1B. Priors for random intercepts landscape
    for (i in 1:NumLands) {a[i] ~ dnorm(0, tau.La) } 

    #1D. Priors for variances for random intercepts
    #    In case of non-mixing, change to half-Cauchy(25) or 
    #    half-Cauchy(16)
    sigma.La ~ dunif(0, 5)
    tau.La  <- 1 / (sigma.La * sigma.La)
    
    #1E. Half-cauchy(25) prior for size
    num.s   ~ dnorm(0, 0.0016)       #<----from NB GLM exercise
    denom.s ~ dnorm(0, 1)            #<----from NB GLM exercise
    size <- abs(num.s / denom.s)     #<----from NB GLM exercise


    ###################################
    #2. Likelihood
    for (i in 1:N) {
       Y[i] ~  dnegbin(p[i], size)
       p[i] <- size / (size + mu[i])  
    	
       #log-link 
       log(mu[i]) <- eta[i]
       eta[i]     <- inprod(beta[], X[i,]) + a[reLands[i]]  

                                          
       #3. Discrepancy measures: 
       #   Expected values, variance, Pearson residuals
       ExpY[i] <- mu[i] 
       VarY[i] <- mu[i] + mu[i] * mu[i] / size
       PRes[i] <- (Y[i] - ExpY[i]) / sqrt(VarY[i])   
       
       #Simulate data from a ZIP  
       YNew[i] ~  dnegbin(p[i], size) 

       #Pearson residuals for simulated data
       PResNew[i] <- (YNew[i] - ExpY[i]) / sqrt(VarY[i]) 
       D[i]       <- pow(PRes[i], 2)                      
       DNew[i]    <- pow(PResNew[i], 2)   
    }     
     Fit         <- sum(D[1:N])
     FitNew      <- sum(DNew[1:N])
}
",fill = TRUE)
sink()
#####################################
#


#Inits function
#Set initial values
inits  <- function () {
  list(beta      = rnorm(K, 0, 0.01), # betas 
       a         = rnorm(NumLands, 0, 0.1), #random effects
       sigma.La  = runif(1, 0, 5),    #sigma landscape
       num.s     = rnorm(1, 0, 25),   #for size
       denom.s   = rnorm(1, 0, 1)    #for size
       )  }


#Step 4: Parameters to estimate
params <- c("beta", "size", 
            "sigma.La", 
            "Fit", "FitNew", "PRes",
            "ExpY"
            )

#Step 5: Run JAGS
NB1   <- jags(data        = win.data,
                inits      = inits,
                parameters = params,
                model      = "NBGLMM.txt",
                n.thin     = 10, 
                n.chains   = 3,
                n.burnin   = 4000,
                n.iter     = 5000)

#And do some updating
NB2  <- update(NB1, n.iter = 10000, n.thin = 10)  
#NB3  <- update(NB2, n.iter = 50000, n.thin = 10)  
#save(NB3, file = "BF_NB3.RData")

#It would be better to do at least another 100K iterations, 
#but this means that not all participants would be able 
#to finish the analysis.

#To reload results from an earlier session, use:
#load(file = "BF_NB3.RData")
#recompile(NB3)
#load.module("dic")

out.NB <- NB2$BUGSoutput
# If you deciced to use NB3, then adjust the code
# above.



#Step 6: Present output 
print(out.NB, digits = 3)  


#Step 7: Assess mixing and check overdispersion
MyNames <- c(colnames(X), "sigma landscape", "size")
#Adjust this variable if extra parameters are added!

MyBUGSChains(out.NB, 
             c(uNames("beta", K), "sigma.La", "size"),
              PanelNames = MyNames)
# We need more iterations here! But that requires half
# an hour computing time


############################################
#Check for overdispersion / underdispersion
#Bayesian approach:
mean(out.NB$sims.list$Fit >  out.NB$sims.list$FitNew) 
#That is perfect

#Plot the sum of squares for each iteration. 
Min <- min(out.NB$sims.list$Fit, out.NB$sims.list$FitNew)
Max <- max(out.NB$sims.list$Fit, out.NB$sims.list$FitNew)

xyplot(out.NB$sims.list$FitNew ~ out.NB$sims.list$Fit, 
       aspect = "iso", 
       col = 1, 
       pch = 16, 
       cex = 0.7,
       xlim = c(Min, Max), 
       ylim = c(Min, Max)) 
#Perfect

#What does the frequentist way of thinking tell us?
Ep <- out.NB$mean$PRes
N <- nrow(BF)
p <- K  + 1 + 1
sum(Ep^2) / (N - p)



#As in the frequentist approach, we now
#need to apply a detailed model validation
#and figure out why we have overdispersion


#Step 8: Numerical output
OUT1 <- MyBUGSOutput(out.NB, 
                     c(uNames("beta", K), "sigma.La", "size"),
                     VarNames = MyNames)
print(OUT1, digits = 5)




#Your homework:
#1. Deal with model selection. We choose these covariates
#   based on some initial analyses.


#Task:
#Apply a model validation


#To understand why we have overdispersion,
#it can also be handy to plot the model fit



#Calculate the fitted values.
#1. Get the betas from the MCMC iterations
#2. Define a grid of covariate values
#3. Calculate the predicted values on this 
#   grid for each MCMC iteration.
#4. Calculate the 95% credible intervals.
#5. Plot the whole thing


#1. Get the betas and gammas
beta.mcmc  <- out.NB$sims.list$beta 
dim(beta.mcmc)

#Using these we can puzzle together the pieces.
#Note that we are not adding the random effects
#to the count part.

#2. Define a grid of covariate values
range(BF$NLowBranchesc)
range(BF$Distance2Edgec)

MyData <- expand.grid(NLowBranchesc  = seq(-2.15, 2.02, length = 25),
                      Distance2Edgec = seq(-0.96, 1.55, length = 25))
                       
#B. Create X matrix with expand.grid
X <- model.matrix(~ 1 + NLowBranchesc + Distance2Edgec, data = MyData)


#C. Calculate the predicted MCMC values
#   for mu
eta.mcmc  <- X %*% t(beta.mcmc) 
mu.mcmc   <- exp(eta.mcmc)



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

ExpY <- MyData2[,"mean"]  #Expected values Y. For a ZIP these will change. 



#So..these are the two variables that we use
#to calculate our grid for:
BD.25   <- seq(-2.15, 2.02, length = 25)
Prc.25  <- seq(-0.96, 1.55, length = 25)
        
#And we convert the vector with expected values into
#a 25 by 25 matrix

ExpY.2d <- matrix(ExpY, nrow = length(Prc.25), ncol = length(BD.25))


#Now plot the jittered data again
library(rgl)
set.seed(12345)
BF$Covariate1.jitter <- jitter(BF$NLowBranchesc, amount = 0.1)
BF$Covariate2.jitter <- jitter(BF$Distance2Edgec, amount = 0.1)

plot3d(x = BF$Covariate1.jitter ,
       y = BF$Covariate2.jitter,
       z = BF$NEggs,
       type = "p",
       size = 10,
       lit = FALSE,
       xlab = "NLowBranchesc",
       ylab = "Distance2Edgec",
       zlab = "Number of eggs"
       )

#Add the surface for the fitted Poisson values
surface3d(BD.25, Prc.25, ExpY.2d, 
          alpha = 0.6, 
          front = "lines", 
          back = "lines", 
          color = "black"
          )



