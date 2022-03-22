#    Introduction to Introduction to MCMC, 
#    Linear Mixed Effects models and GLMM with R
#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
#Spider demo code


######################################################
#set the working directory on a Mac & read the data

Spiders2 <- read.table(file = "Spiders.txt", 
                       header = TRUE, 
                       dec = ".")

names(Spiders2)
str(Spiders2)

###################################################################
#Load packages and library files
library(lattice)  #Needed for multi-panel graphs
#source the file: HighstatLibV6.R
source(file="/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV6.R")  
##################################################################



##################################################################
#Housekeeping
Spiders2$fPlot <- factor(Spiders2$Plot)

#We will drop some sites.
#This will allow us to demonstrate random intercept AND slope
#models later on 
Spiders <- Spiders2[
                    Spiders2$fPlot != "4" &
                    Spiders2$fPlot != "9" &
                    Spiders2$fPlot != "11" &
                    Spiders2$fPlot != "14" &
                    Spiders2$fPlot != "23",]
Spiders$fPlot <- as.factor(as.numeric(Spiders$fPlot)) 
##################################################################



##################################################################
#Data exploration
plot(x = Spiders$HerbLayer, 
     y = Spiders$Hlog10,
     xlab = "Percentage of herb layer",
     ylab = "Shannon index",
     cex.lab = 1.5, 
     pch = 16)

#And add a straight line
M0 <- lm(Hlog10 ~ HerbLayer, 
            data = Spiders)
summary(M0)


#Model validation
#Plot residuals versus Plot
E0 <- resid(M0)
boxplot(E0 ~ fPlot, 
        data = Spiders,
        xlab = "Plot",
        ylab = "Residuals",
        cex.lab = 1.5)
abline(h = 0, lty = 2)

#Return to Powerpoint  <------
######################################################






######################################################
#Block 2
#Add fPlot as covariate
M1 <- lm(Hlog10 ~ HerbLayer + fPlot, 
            data = Spiders)
summary(M1)
drop1(M1, test = "F")


#Plot fitted values
plot(x = Spiders$HerbLayer, 
     y = Spiders$Hlog10,
     xlab = "Percentage of herb layer",
     ylab = "Shannon index",
     cex.lab = 1.5, 
     pch = 16)


PlotLevels <- unique(levels(Spiders$fPlot))
for (i in PlotLevels){
   MinMaxi <- range(Spiders$HerbLayer[Spiders$fPlot==i])
   MyDatai <- data.frame(HerbLayer=
                  seq(from = MinMaxi[1],
                      to   = MinMaxi[2],
                      length=10),
                  fPlot = i)
   Pi <- predict(M1, newdata = MyDatai)
   lines(MyDatai$HerbLayer, Pi, lty=1,col=1,lwd=1)
}

#Return to Powerpoint <---------
#############################################################





#############################################################
#Block 3
library(lme4)    #You need to download this
#fPlot needs to be a factor

M2 <- lmer(Hlog10 ~ HerbLayer + (1 | fPlot), 
            data = Spiders)
summary(M2)

icc <- 0.1788/(1.788 + 0.01524)

#REML: ‘restricted maximum likelihood.’ 
#It is a modified version of maximum likelihood estimation 
#and provides better estimates for the variance terms. 
#Further information on ML and REML can be found in Zuur et al. (2009a).

#Estimated variances and square roots of variances
#Fixed effects

#No p-values
#Bates: 
# t-statistic in lmm does not follow a t-distribution
# F-statistic in lmm does not follow an F-distribution
#...whatever the degrees of freedom.
#Especially not for small samples

#You want to have p-values? Like
#other packages?
Betas  <- fixef(M2)               #Get the betas
SE     <- sqrt(diag(vcov(M2)))    #Get the SEs
tVal   <- Betas  / SE
pval   <- 2*pnorm(-abs(tVal))     #Z distribution
Output <- cbind(Betas, SE, tVal, pval)
print(Output, digits = 3)
#Voila...there are your p-values
#Use them with care!

################################################
#To understand what this model is doing:
#Write down the equation and sketch fitted values

#Equation:
summary(M2)
#Shannon_ij = 0.91 + 0.16 * HerbLayer_ij + a_i + eps_ij

#Fixed part:
#    E(Shannon_ij) = 0.91 + 0.16 * HerbLayer_ij

#Random part:
#     a_i + eps_ij


#Let us sketch the fixed part
#1. Create covariates on a grid
#2. Get X matrix for these covariates
#3. Extract betas from M2
#4. Calculate mu = X * beta
#5. Plot mu versus HerbLayer
#6. Add CIs....

#1. Creat covariates on a grid
MinMax <- range(Spiders$HerbLayer)
MyData <- data.frame(HerbLayer=
                      seq(from   = MinMax[1],
                          to     = MinMax[2],
                          length = 10))

#2. Convert these into a matrix X
X <- model.matrix(~ HerbLayer, data = MyData)

#3. Extract betas
betas <- fixef(M2)

#4. Calculate the fitted values
mu <- X %*% fixef(M2)

#5. Plot mu versus HerbLayer
#We will do 6. in the exercises

plot(x = Spiders$HerbLayer, 
     y = Spiders$Hlog10,
     xlab = "Percentage of herb layer",
     ylab = "Shannon index",
     cex.lab = 1.5, pch = 1)
lines(MyData$HerbLayer, mu, lwd = 5)

#That is the fixed part of the model.
#And for a paper only add the CIs, and you are done

#For understanding what the mm is doing, 
#we show the effect of the random effects a_i
#For each plot we need to add an a_i to the line
#Let's get the random intercepts first
ranef(M2)
a <- ranef(M2)$fPlot$'(Intercept)'
a
#If you have a lot of them, then check them for
#normality and independence (e.g. in space)

#This is the data and the fixed part of the model
#again. Now add the random intercepts.
plot(x = Spiders$HerbLayer, 
     y = Spiders$Hlog10,
     xlab = "Percentage of herb layer",
     ylab = "Shannon index",
     cex.lab = 1.5, pch = 1)

lines(MyData$HerbLayer, mu, lwd = 5, col = 2)

#Use a simple loop to add the a_i effect to the
#fitted values. But....the range of HerbLayer cover 
#differs per plot....hence the reason the loop looks
#not so nice

j<-1
PlotLevels <- unique(Spiders$fPlot)
for (i in PlotLevels){
   MinMaxi <- range(Spiders$HerbLayer[Spiders$fPlot==i])
   MyDatai <- data.frame(HerbLayer=
                  seq(from   = MinMaxi[1],
                      to     = MinMaxi[2],
                      length = 10),
                  fPlot = i)
   X <- model.matrix(~ HerbLayer, data = MyDatai)
               
   Pi <- X %*% betas + a[j]
   j <- j + 1
   lines(MyDatai$HerbLayer, Pi, lty=1,col=1,lwd=1)
}

#Question 1 for you:
#Where are the epsilons?



#Question 2 for you.
#This is the model
#Shannon_ij = 0.91 + 0.16 * HerbLayer_ij + a_i + eps_ij
#a_i ~ N(0, sigma^2_Plot)

#What is the effect of sigma^2_Plot ??
#Formulated differently:
#   Suppose that sigma^2_Plot is small.
#   How do the fitted values look like?
#And what is  sigma^2_Plot is big?

#Simulate some examples.
#Run the code with each of these:
#The estimated sigma_Plot is 0.1337
a.simulated <- rnorm(26, mean = 0, sd = 0.01)
#a.simulated <- rnorm(26, mean = 0, sd = 0.5)

plot(x = Spiders$HerbLayer, 
     y = Spiders$Hlog10,
     xlab = "Percentage of herb layer",
     ylab = "Shannon index",
     cex.lab = 1.5, pch = 1)

lines(MyData$HerbLayer, mu, lwd = 5, col = 2)

j<-1
PlotLevels <- unique(Spiders$fPlot)
for (i in PlotLevels){
   MinMaxi <- range(Spiders$HerbLayer[Spiders$fPlot==i])
   MyDatai <- data.frame(HerbLayer=
                  seq(from   = MinMaxi[1],
                      to     = MinMaxi[2],
                      length = 10),
                  fPlot = i)
   X <- model.matrix(~ HerbLayer, data = MyDatai)
               
   Pi <- X %*% betas + a.simulated[j]
   j <- j + 1
   lines(MyDatai$HerbLayer, Pi, lty=1,col=1,lwd=1)
}
#################################################
#Return to Powerpoint







#############################################################
#Block 4


#In order to compare results with JAGS later,
#standardize covariates

#Standardizing may also help for GLMMs fitted
#with lme4!


#Standardize covariates:
MyNorm <- function(x){ (x - mean(x)) / sd(x)}
Spiders$HerbLayerc <- MyNorm(Spiders$HerbLayer)
Spiders$GroundVegc <- MyNorm(Spiders$GroundVeg)
Spiders$Litterc    <- MyNorm(Spiders$Litter)

#Step 1: By design of the study use plot as random effect
#And use the three continuous covariates
M3 <- lmer(Hlog10 ~ HerbLayerc + GroundVegc + Litterc + 
           (1 | fPlot), 
           data = Spiders)

#Step 2: Are all covariates significant?
summary(M3) #Get p-values for the t-values? 
            #See code above
            #This is done with REML
            #We called it 'step 7A' in the old protocol


# Or apply a likelihood ratio test. 
# We called it 'step 7c' in the old protocol
# But this likelihood ratio test needs ML estimation..not REML
# The 'Why' question can only be answered with difficult 
# mathematics
M4 <- lmer(Hlog10 ~ HerbLayerc + GroundVegc + Litterc + 
           (1 | fPlot), 
            data = Spiders, REML = FALSE)

M4A <- update(M4, .~. - HerbLayerc)
M4B <- update(M4, .~. - GroundVegc)
M4C <- update(M4, .~. - Litterc)
anova(M4, M4A)
anova(M4, M4B)
anova(M4, M4C)

#This is exactly the same as:
drop1(M4, test = "Chi")
#So...better use this one (but drop1 does not work with 
#                          the package nlme!)

#Or use AIC:
AIC(M4, M4A, M4B, M4C)
#But use AIC with ML!!!!
#And don't mix REML and ML AICs!!

#Either drop GroundVegc...or keep it in
#That depends on your own attitude to model selection

#Step 3: Present the final model..and explain what it all means
#Use REML for this!!
M5 <- lmer(Hlog10 ~ HerbLayerc + GroundVegc + Litterc + 
           (1 | fPlot), 
           data = Spiders, REML = TRUE)

summary(M5)


#Check whether this model does not violate anything!
#Homogeneity
#Independence
#No non-linear patterns in the residuas
#Normality of residuals

#Get residuals and fitted values
E5 <- resid(M5)  
F5 <- fitted(M5)
#In the exercise we will see what exactly this is

#Check homogeneity
plot(x = F5, 
     y = E5,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)
#That is ok

#Make variograms of the residuals
#But that is a different course.

#Check for non-linear patterns in the residuals
plot(x = Spiders$HerbLayerc, 
     y = E5,
     xlab = "HerbLayerc",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)
#That is ok...and do it for all (!) other available
#covariates.

#Check for normality
hist(E5)


#Try to explain what it all means
#Write down the fitted model
summary(M5)


#And sketch the fitted values.
# 1. Specify covariate values for which to do predictions.
# 2. Convert the covariate values to a matrix X using the 
#    expand.grid function. Extract the estimated regression 
#    parameters with the fixef function.
# 3. Calculate the predicted values via X × β.  
# 4. Calculate standard errors for the predicted values by 
#    taking the square root of the diagonal elements of 
#    X × cov(β) × X’.
# 5. Plot predicted values versus a covariate.
# 6. Add 1.96 ± the standard errors to the fitted values.



range(Spiders$HerbLayerc)
MyData <- data.frame(HerbLayerc = seq(-1.3, 2, length = 10),
                     GroundVegc = 0,  #= average 
                     Litterc = 0)     #= average
X         <- model.matrix(~HerbLayerc + GroundVegc + Litterc,  
                          data = MyData)
Betas     <- fixef(M5)
FitManual <- X %*% Betas
SE        <- sqrt(diag(X %*% vcov(M5) %*% t(X)))

plot(x = MyData$HerbLayerc,
     y = FitManual,
     type = "l",
     xlab = "Herb layer cover (standardized)",
     ylab = "Fitted values for average GrounVeg and Litter values")
lines(x = MyData$HerbLayerc,
      y = FitManual + 1.96 * SE,
      lty = 2)
lines(x = MyData$HerbLayerc,
      y = FitManual - 1.96 * SE,
      lty = 2)
      
           
#For an R2, see:
#Nakagawa and Schielzeth (2013)
#A general and simple method for obtaining R2 from generalized linear mixed-effects models
#Methods in Ecology and Evolution
#Volume 4, Issue 2, pages 133–142, February 2013
######################################################

#Return to powerpoint	





####################################################
#Block 5
#Run a mixed model in JAGS
#Load the R2jags package and the support routine 
library(R2jags)
source("/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstat.R")


#Prepare data
X <- model.matrix(~ HerbLayerc + GroundVegc + Litterc, 
                    data = Spiders)
K <- ncol(X)
head(X)
K

#We need a vector numbers for each plot
#Needs to start at 1
Plot <- as.numeric(as.factor(Spiders$fPlot))
Nre  <- length(unique(Spiders$fPlot))
Plot
Nre  # Number of random intercepts

#Get all data
win.data <- list(Y    = Spiders$Hlog10,
                 X    = X,
                 K    = K,
                 N    = nrow(Spiders),
                 Plot = Plot,
                 Nre  = Nre)
win.data
#Check the flowchart in Powerpoint to see where we are



#Specify modelling code for JAGS
#Modelling code
sink("lmm.txt")
cat("
model {
    #1A. Priors for regression parameters
    for (i in 1:K) {beta[i] ~ dnorm(0, 0.0001)} 

    #1B. Priors for random effect plot
    for (i in 1:Nre) {a[i] ~ dnorm(0,  tau.plot)}
    
    #1C. Priors for sigma_plot
    tau.plot <- 1 / (sigma.plot * sigma.plot)
    sigma.plot ~ dunif(0.001, 5)
    
    #1D. Priors for sigma_eps
    tau.eps <- 1 / (sigma.eps * sigma.eps)
    sigma.eps ~ dunif(0.001, 5)
    
    #2. Likelihood
    for (i in 1:N) {
        Y[i]  ~  dnorm(mu[i], tau.eps)
        mu[i]  <- eta[i] 
        eta[i] <- inprod(beta[], X[i,]) + a[Plot[i]]
      }
    }
",fill = TRUE)
sink()



#Initial values for each chain, for each parameter
inits <- function () {
   list(beta       = rnorm(K, 0, 0.01),
        a          = rnorm(Nre, 0, 1),
        sigma.eps  = runif(1, 0.001, 5),
        sigma.plot = runif(1, 0.001, 5)
        )}

#Parameters to sore
params <- c("beta", "a", 
            "sigma.plot", 
            "sigma.eps")

#Run JAGS
J0 <- jags(data       = win.data,
           inits      = inits,
           parameters = params,
           model.file = "lmm.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)
             
J1  <- update(J0, n.iter = 10000)
out <- J1$BUGSoutput
print(out, digits = 4)

#Assess mixing
MyBUGSChains(out, c(uNames("beta", K),
                    "sigma.plot", 
                    "sigma.eps"))

#Present output
OUT1 <- MyBUGSOutput(out, c(uNames("beta", K), 
                          "sigma.plot", 
                          "sigma.eps"))
print(OUT1, digits =3)
#Compare with lme4
summary(M5)
#That is similar.


#Visual presentation
library(coefplot2) 
coefplot2(OUT1[1:4,1], 0*OUT1[1:4,2], offset = 0, col = 1,
          varnames = colnames(X), cex.var = 1,
          xlim = c(-0.1,1),
          lower2 = OUT1[1:4,3],
          upper2 = OUT1[1:4,4], main ="MCMC results")

#M5 results
#Get betas and ses from M5:
beta5 <- fixef(M5)
se5   <- sqrt(diag(vcov(M5)))
coefplot2(beta5, 0*se5, offset = 0.05, col = 1,
          varnames = names(beta5), 
          cex.var = 1,
          xlim = c(-0.1,1),
          lower2 = beta5 - 1.96*se5,
          upper2 = beta5 + 1.96*se5, add= TRUE)





         
         
         
         
         
         


#In the exercises we will:
#1. Show how to do model selection with MCMC
#2. Show how to do model validation with MCMC
#3. Show how to calculate fitted values, residuals, etc.
#4. Present model results.

#You will also see 12 more examples of JAGS...so plenty of 
#opportunities to fully understand the code
#####################################################
#Return to powerpoint



