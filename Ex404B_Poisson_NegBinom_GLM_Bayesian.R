#    Highland Statistics Ltd.
#    www.highstat.com
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  



######################################################
#Squirrel example
#Set the working directory and import the data

setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data")
SQ <- read.table(file   = "RedSquirrelsV2.txt", 
                 header = TRUE, 
                 dec    = ".")

str(SQ)
names(SQ)
###################################################################




###################################################################
#Load packages and library files
library(lattice)  #Needed for multi-panel graphs
library(R2jags)
source(file = "HighstatLibV9.R")  
source(file = "MCMCSupportHighstatV4.R")
##################################################################


##################################################################
#Data exploration
#Outliers
MyVar <- c("SqCones", "Ntrees", "DBH", 
           "TreeHeight", "CanopyCover")

Mydotplot(SQ[,MyVar])
#Ouch....1 large DBH value!
#A few small Canopy values!
#A few large Ntrees values.
#  There may be a patch with large
#  DBH and NTrees value!

#Remove the value with extreme large DBH value
SQ2 <- subset(SQ, DBH < 0.6)
dim(SQ)
dim(SQ2)


#Collinearity
MyVar <- c("Ntrees", "DBH", 
           "TreeHeight", 
           "CanopyCover")

pairs(SQ2[, MyVar], 
      lower.panel = panel.cor)

corvif(SQ2[,MyVar])
#In principle that is ok!!
#But the outliers may be doing funny things!


#Relationships
MyVar <- c("Ntrees", "DBH", 
           "TreeHeight", 
           "CanopyCover")
Myxyplot(SQ2, MyVar, "SqCones", 
         MyYlab = "Number of stripped cones")
#Yes...that point was indeed going to cause trouble!!!


#Zero inflation?
plot(table(SQ2$SqCones))
#No

#Conclusion data exploration: 
#1 outlier in DBH
#Small amount of collinearity
#Weak patterns between Y and X
#Potentially NB distribution for SqCones
#####################################################





#####################################################
#Start Bayesian analysis

#Standardize continuous covariates
MyStd <- function(x) { (x - mean(x)) / sd(x)}

SQ2$Ntrees.std      <- MyStd(SQ2$Ntrees)
SQ2$TreeHeight.std  <- MyStd(SQ2$TreeHeight)
SQ2$CanopyCover.std <- MyStd(SQ2$CanopyCover)


######################################
#This is the reference model
M1 <- glm(SqCones ~ Ntrees.std + TreeHeight.std +
                    CanopyCover.std,
          family = "poisson",
          data = SQ2)


##################################################                     
#Use MCMC to estimate the parameters                     
#1. Bundle data

X <- model.matrix(~ Ntrees.std +  
                    TreeHeight.std + 
                    CanopyCover.std, 
                  data = SQ2)
   
                   
K <- ncol(X)  #Number of columns
head(X)


win.data <- list(Y    = SQ2$SqCones,  #Response 
                 X    = X,            #Covariates
                 N    = nrow(SQ2),    #Sample size
                 K    = K)            #Number of betas
win.data

###################################################
# 2. JAGS modelling code
# In each JAGS exercise we will add something new.
# In this exercise we show two new things
#  A. How to do a Poisson GLM in JAGS
#  B. How to assess overdispersion in MCMC
#  C. How to do a negative binomial GLM in JAGS

#We will use a model with CanopyCover. Model selection
#is discussed later.

#We will start with A



############################################################
#This is the JAGS code for a Poisson GLM
#It is actually much easier than that of a linear regression model!
sink("squirrelsGLMJAGS.txt")
cat("
model{
    #1A. Priors beta 
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}

    #2. Likelihood
    for (i in 1:N) {
      Y[i]        ~ dpois(mu[i])      
      log(mu[i]) <- inprod(beta[], X[i,]) 
      # or <- beta[1] * X[i, 1] + beta[2] * X[i, 2] + beta[3] * X[i, 3] + beta[4] * X[i, 4]
     
      #3. Discrepancy measures 
      #Pearson residuals (normalized)

      Exp[i] <- mu[i] 
      Var[i] <- mu[i]
      E[i]   <- (Y[i]  - Exp[i]) / sqrt(Var[i])    
  
      #Simulated data with mean/variance taken from the fitted model 
      #See text under block B, below.  

      YNew[i] ~  dpois(mu[i])                      

      #Pearson residual for predicted data      
      ENew[i] <- (YNew[i] - Exp[i]) / sqrt(Var[i]) 
      
      #Squared residuals
      D[i]    <- E[i]^2                      
      DNew[i] <- ENew[i]^2   
     }          

     #Sum of squared Pearson residuals:
     Fit     <- sum(D[1:N])      
     #Sum of squared predicted Pearson residuals:     
     FitNew  <- sum(DNew[1:N])   
}
",fill = TRUE)
sink()
#####################################
#

 
#3. Initial values & parameters to save
inits  <- function () {
  list(
    beta = rnorm(K, 0, 0.01))  }

params <- c("beta", "E", "mu",  
            "Fit", "FitNew")


#4. Start JAGS
J1 <- jags(data       = win.data,
           inits      = inits,
           parameters = params,
           model      = "squirrelsGLMJAGS.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

J2  <- update(J1, n.iter = 10000, n.thin = 10)  
out <- J2$BUGSoutput
print(J2, digits = 3)  



#5. Assess mixing
MyNames <- c("Intercept", "Number of Trees", 
             "Tree height", "Canopy cover")


MyBUGSChains(out, c(uNames("beta", K)), PanelNames = MyNames)
MyBUGSACF(out, c(uNames("beta", K)), PanelNames = MyNames)

#6. Present output
OUT1 <- MyBUGSOutput(out, c(uNames("beta", K)),
                     VarNames = MyNames)
print(OUT1, digits =5)
MyBUGSHist(out, c(uNames("beta", K)), PanelNames = MyNames)


OUT1
summary(M1)  #Small differences...but for such a small
             #data set the glm results may also be off!        
######################################################


##################################################
#B: Check for overdispersion (Frequentist way of thinking)

#Option 1:
E    <- out$mean$E
N    <- nrow(SQ2)
p    <- K
sum(E^2) / (N - p)


#Option 2:
# In each JAGS exercise we will add something new.
# In this exercise we add tools to assess the quality of the fit.
# Recall that in JAGS we create a large number of iterations
# of the parameters. In the code above we simulated data
# from the Poisson model. If the simulated data is similar 
# to the observed data then we have a good model. One way 
# to quantify "similar" is by taking:
#   the sum of squared Pearson residuals for the observed data
#   the sum of squared Pearson residuals for the simulated data
# and compare these.... (other options exist, see Lawson (2013)).
# Using the large number of realisations of these two statistics, we
# can count how often the residual sum of squares for the simulated 
# data is larger than that for the original data. Hopefully it 
# is 50-50. That would mean that the simulated data is indeed 
# similar to the observed data. This leads to a Bayesian p-value.

mean(out$sims.list$Fit > out$sims.list$FitNew)

#Fitted values always have larger residual sum of squares
#That means overdispersion

# You can also evaluate the occurance of zeros, maximum values



#You now need to assess why there is overdispersion
#See the frequentist version of this exercise:
# A. outliers
# B. etc

#Plot E vs each covariate
#Plot E vs fitted values
#We leave this as an exercise
##########################################################






#####################################################################
#Apply a NB GLM  
#We now show code for fitting a NB GLM in JAGS.
#We are going to use the optimal model (saving you
#some time on model selection)


#1. Bundle data
X <- model.matrix(~ CanopyCover.std, data = SQ2)
K <- ncol(X)  #Number of columns
head(X)

win.data <- list(Y    = SQ2$SqCones, #Response
                 X    = X,           #Intercept + covariate
                 N    = nrow(SQ2),   #Sample size
                 K    = K)           #Number of betas
win.data

###################################################
# 2. JAGS modelling code

# Y_i ~ negbinom(mu_i, size)
# E(Y_i) = mu_i = exp(model)
# var(Y_i) = mu_i + mu_i^2 / size


sink("squirrels_NB_GLMJAGS.txt")
cat("
model{
    # 1A. Priors for beta
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}

    # 1B. Uniform prior for size (= k = theta = 1 /alpha = alpha)
    # restrict size parameter to low number, and above zero (same for sigmas)
    size ~ dunif(0.0001, 5)  #k or theta (second neg binom parameter)

    # 2. Likelihood
    for (i in 1:N) {

      # This is just how JAGS implements a NB distribution
      # It is based on general mathematical rules for a NB

      Y[i] ~  dnegbin(p[i], size)
      p[i] <- size / (size + mu[i])  
      
      log(mu[i]) <- eta[i]  
      eta[i]     <- inprod(beta[], X[i,]) 
                    #= beta[1] * X[i,1] + .. + beta[K] * X[i,K]
  
      #3. Discrepancy measures 
      #Pearson residuals
      Exp[i] <- mu[i] 
      Var[i] <- mu[i] + mu[i] * mu[i] / size
      E[i]   <- (Y[i]  - Exp[i]) / sqrt(Var[i])    

      #Simulated data with mean/variance taken from the fitted model     
      YNew[i] ~  dnegbin(p[i], size)                      

      #Pearson residual for predicted data      
      ENew[i] <- (YNew[i] - Exp[i]) / sqrt(Var[i]) 
      
      #Squared residuals
      D[i]    <- E[i]^2                      
      DNew[i] <- ENew[i]^2   
     }          
     
     #Sum of squared Pearson residuals and
     #squared predicted Pearson residuals
     Fit         <- sum(D[1:N])        
     FitNew      <- sum(DNew[1:N])   
}
",fill = TRUE)
sink()
#####################################

#3. Initial values & parameters to save
inits  <- function () {
  list(
    beta     = rnorm(K, 0, 0.01),
    size     = runif(1,0,5 ) 
    )  }

params <- c("beta", "E", "mu", "size",  
            "Fit", "FitNew")


#4. Start JAGS
NB1 <- jags(data       = win.data,
           inits      = inits,
           parameters = params,
           model      = "squirrels_NB_GLMJAGS.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

NB2    <- update(NB1, n.iter = 10000, n.thin = 10)  
out.nb <- NB2$BUGSoutput
print(NB2, digits = 3)  



#5. Assess mixing
MyNames <- c("Intercept", "Canopy cover", "size")


MyBUGSChains(out.nb, c(uNames("beta", K), "size"), 
             PanelNames = MyNames)
MyBUGSACF(out.nb, c(uNames("beta", K), "size"), 
          PanelNames = MyNames)

#6. Present output
OUT1 <- MyBUGSOutput(out.nb, c(uNames("beta", K), "size"), 
                     VarNames = MyNames)
print(OUT1, digits =5)
MyBUGSHist(out.nb, c(uNames("beta", K), "size"), PanelNames = MyNames)


OUT1
#Compare with:
library(MASS)
M2 <- glm.nb(SqCones ~ CanopyCover.std,
             data = SQ2)
summary(M2)
#That is good
######################################################

#Check for overdispersion
#Option 1:
E    <- out.nb$mean$E
N    <- nrow(SQ2)
p    <- K + 1
sum(E^2) / (N - p)


#Option 2:
mean(out.nb$sims.list$Fit > out.nb$sims.list$FitNew)
#That is all good!

#Model validation is similar as before




######################################################
#Model presentation
#Time for some ggplot2 coding!

library(ggplot2)


range(SQ2$CanopyCover.std)
MyData <- data.frame(CanopyCover.std = seq(from = -3.95 , to = 0.89 , 
                                            length = 25))

Xp     <- model.matrix(~ CanopyCover.std, data = MyData)


#Get the MCMC betas
Beta.mcmc <- out.nb$sims.list$beta 
Beta.mcmc

#Calculate the fitted values
eta.mcmc <- Xp %*% t(Beta.mcmc)
mu.mcmc  <- exp(eta.mcmc)  
  
#Note that this mu.mcmc contains
#predictions for the 25 artifical covariate values


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
head(MyData2)



# Each row in L is an artificial depth value
# The first column is the posterior mean of the MCMC iterations
# at that specific fitted value.
# The second row is the SE, the third and fourth are the 
# 2.5 and 97.5% quartiles and are used for the 95% credible intervale
# Now plot the data and draw lines for the mean,
# and 2.5% and 97.5% values:


#Calculate the covariate on the original scale
MyData2$CanopyCover <- MyData2$CanopyCover.std * sd(SQ2$CanopyCover) + 
                       mean(SQ2$CanopyCover)
#This is a simple back transformation

p <- ggplot()
p <- p + geom_point(data = SQ2, 
                    aes(y = SqCones, x = CanopyCover),
                    shape = 16, 
                    size = 1.5)
p <- p + xlab("CanopyCover") + ylab("Number of stripped cones")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData2, 
                   aes(x = CanopyCover, 
                       y = mean))

p <- p + geom_ribbon(data = MyData2, 
                     aes(x = CanopyCover, 
                         ymax = up, 
                         ymin = lo),
                         alpha = 0.2)
p
###########################################################
