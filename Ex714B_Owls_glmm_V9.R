#    Highland Statistics Ltd.
#    www.highstat.com
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


####################################################
#Import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/FilesMCMC_GLMM/CDMCMCGLMGLMMCourse/Data")
Owls <- read.table("Owls.txt", header=TRUE)


#Check imported data
names(Owls)
str(Owls)

#This is what we have
# [1] "Nest"               "Xcoord"             "Ycoord"            
# [4] "FoodTreatment"      "SexParent"          "ArrivalTime"       
# [7] "SiblingNegotiation" "BroodSize"          "NegPerChick"       
# [10] "NCalls"          


#Task: Model SiblingNegotiation  as a function of:
#  -FoodTreatment    (factor)   
#  -SexParent         (factor)
#  -ArrivalTime       (continuous)
#  -FoodTreatment x SexParent
#  -SexParent x ArrivalTime
#  -Nest effect (but nests were selected randomly from a large group)

###################################################################
#Load packages and library files
library(lattice)  #Needed for multi-panel graphs

#source the file: HighstatLibV9.R
source(file="/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV9.R")  
#It can be dowloaded from the course website. During this course we will extend it.

library(R2jags)
source("/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstatV4.R")

##################################################################


#SiblingNegotiation is too long....use shorter name:
Owls$NCalls <- Owls$SiblingNegotiation      



##########################################
#Data exploration
boxplot(NCalls ~ Nest, data = Owls, varwidth = TRUE)
#Nest effect?


#How many nests, and how many observations per nest do we have?
length(unique(Owls$Nest))

table(Owls$Nest)
plot(table(Owls$Nest), type = "h")


# #Potential trouble with 2 nests with only 4 observations
# #Remove these
# Owls2 <- Owls[Owls$Nest!="Chevroux" & Owls$Nest!="Forel" & 
              # Owls$Nest!="GDLV" & Owls$Nest!="SEvaz",]
# Owls2$Nest <- factor(Owls2$Nest)

Owls2 <- Owls

#What is the percentage of zeros?
100 * sum(Owls$NCalls == 0, na.rm = TRUE) / nrow(Owls)
plot(table(Owls$NCalls), type = "h")


#Look at some design and interaction plots
plot.design(NCalls ~ FoodTreatment + SexParent + Nest ,
            data = Owls2)
            
interaction.plot(Owls2$FoodTreatment, Owls2$SexParent,
                 Owls2$NCalls, xlab="FoodTreatment",
                 ylab="Sibbling negotiation")


#Plot the data as time series
library(lattice)
xyplot(NCalls ~ ArrivalTime | Nest, type="h",
       data = Owls2, col = 1)
  
xyplot(NCalls ~ ArrivalTime | Nest, type="h",
         subset = (FoodTreatment == "Satiated"),
         data = Owls2, col = 1)

xyplot(NCalls ~ ArrivalTime | Nest, type="h",
       subset = (FoodTreatment == "Deprived"),
       data = Owls2, col = 1)



xyplot(NCalls ~ ArrivalTime | SexParent * FoodTreatment,
       data = Owls2,
       ylab = "Sibbling negotiation", xlab = "Arrival time (hours)",
         panel=function(x, y){
           panel.grid(h = -1, v = 2)
           panel.points(x, y, col = 1)
           panel.loess(x, y, span = 0.5, col = 1, lwd = 2)})


#Plot the spatial position of the sites
xyplot(Ycoord ~ Xcoord, aspect = "iso", col = 1,
       xlab = "X coordinate",
       ylab = "Y coordinate",
       data = Owls2,
       pch = 16)
#################################################
#For the offset we need:
Owls$LogBroodSize <- log(Owls$BroodSize)


# Fit the following model:
# NCalls_ij ~ Poisson(mu_ij)
# mu_ij  = exp(eta_ij)
# eta_ij = Intercept + SexParent + FoodTreatment + ArrivalTime +
         # SexParent x FoodTreatment +
         # SexParent x ArrivalTime + 1 * LogBroodSize +
         # a_i
# a_i ~ N(0, sigma_Nest^2)          
# i is nest index




########################################################################################
#MCMC part of the exercise


#This will be beneficial in JAGS:
Owls$ArrivalTime.std <- (Owls$ArrivalTime - mean(Owls$ArrivalTime)) / 
                         sd(Owls$ArrivalTime)


#You can either start from scratch with the full model...run different models and 
#use DIC to do model validation, or we can save ourselves some time...pretend that
#we have done all the sub-models...and continue/start with the optimal model:

#Use this for starting values
M6 <- glm(NCalls ~ FoodTreatment + ArrivalTime.std +
                   offset(LogBroodSize),
            family = poisson, 
            data = Owls)





####################################
##Step 1: Fit the model.
#Use MCMC to estimate the parameters                     
#1. Bundle data
X <- model.matrix(~ FoodTreatment + ArrivalTime.std, data = Owls)                      
K <- ncol(X)  #Number of columns
head(X)

#Random effects:
Nest <- as.numeric(as.factor(Owls$Nest))
Nest
Nre <- length(unique(Owls$Nest))


#The code below is copied an pasted from the owls mixed modelling code
#Nest was changed into Hive
win.data <- list(Y     = Owls$NCalls, 
                 X     = X,
                 N     = nrow(Owls),
                 K     = K,
                 Nest  = Nest,
                 Nre   = Nre,
                 LogBS = Owls$LogBroodSize)
win.data

###################################################
# 2. JAGS modelling code

#This code was copied from the coral reef example. We only changed Site into Nest.
sink("OwlsGLMM.txt")
cat("
model{
    #1A. Priors beta and sigma
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}

    #1B. Priors random effects and sigma_Nest
    for (i in 1:Nre) { a[i] ~ dnorm(0, tau_Nest)}
    
    tau_Nest   <- 1 / (sigma_Nest * sigma_Nest)
    num         ~ dnorm(0, 0.0016)       #<----half-Cauchy(25)
    denom       ~ dnorm(0, 1)            #<----half-Cauchy(25)
    sigma_Nest <- abs(num / denom)       #<----half-Cauchy(25)


    #2. Likelihood
    for (i in 1:N) {
      Y[i]        ~ dpois(mu[i])  
      log(mu[i]) <- eta[i] 
      eta[i]     <- inprod(beta[], X[i,]) + a[Nest[i]] + LogBS[i] 
  
      #3. New stuff: Discrepancy measures 
      Exp[i] <- mu[i] 
      Var[i] <- mu[i]
      
      #Pearson residuals
      E[i]   <- (Y[i]  - Exp[i]) / sqrt(Var[i])   
       
      #Simulated Poisson data 
      YNew[i] ~  dpois(mu[i]) 
      
      #Pearson residuals for simulated data
      ENew[i] <- (YNew[i] - Exp[i]) / sqrt(Var[i]) 
      D[i]    <- pow(E[i], 2)                      
      DNew[i] <- pow(ENew[i], 2)   
     }          
     Fit         <- sum(D[1:N])     #Sum of squared residuals  
     FitNew      <- sum(DNew[1:N])  #Sum of squared predicted residuals
}
",fill = TRUE)
sink()
#####################################

#3. Initial values & parameters to save
inits  <- function () {
  list(
    beta  = rnorm(K, coef(M6), 0.1),
    a     = rnorm(Nre, 0, 0.1),
    num   = rnorm(1, 0, 25),   #<-- half-Cauchy(25)
    denom = rnorm(1, 0, 1)  )  #<-- half-Cauchy(25)
    }

params <- c("beta", "E", "a", 
            "sigma_Nest", "mu", 
            "Fit", "FitNew")


#Don't forget to change the file name!
#4. Start JAGS
P1 <- jags(data       = win.data,
           inits      = inits,
           parameters = params,
           model      = "OwlsGLMM.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

P2  <- update(P1, n.iter = 10000, n.thin = 10)  
outP <- P2$BUGSoutput

print(P2, digits = 3)  


#5. Assess mixing
MyBUGSChains(outP, c(uNames("beta", K), "sigma_Nest"))
MyBUGSACF(outP, c(uNames("beta", K), "sigma_Nest"))
#

#6. Present output
OUTP <- MyBUGSOutput(outP, c(uNames("beta", K), "sigma_Nest"))
print(OUTP, digits =5)
MyBUGSHist(outP, c(uNames("beta", K), "sigma_Nest"))
######################################################



######################################################
# Assess goodness of fit: Bayesion p-value:
# We have a large number of FitNew and Fit values. 
# Count how often one is bigger than the other:

mean(outP$sims.list$Fit >  outP$sims.list$FitNew)
#Yes...we expected that based on the frequentist results!

E   <- outP$mean$E
sum(E^2) / (nrow(Owls) - 4)


#Plot the sum of squares for each iteration. 
Min <- min(outP$sims.list$Fit, outP$sims.list$FitNew)
Max <- max(outP$sims.list$Fit, outP$sims.list$FitNew)

xyplot(outP$sims.list$FitNew ~ outP$sims.list$Fit, 
       aspect = "iso", 
       col = 1, 
       pch = 16, 
       cex = 0.7,
       xlim = c(Min, Max), 
       ylim = c(Min, Max)) 

#################################################################
#Same story...why do we have overdispersion?

#Copy paste from the frequentist analysis:
#Why do we have overdispersion?
  #A. Outliers?                  ==> Remove them?
  #B. Missing covariates or interactions?  ==> Go back or..add a latent variable 
  #C. Zero inflation?            ==> ZIP
  #D. Large variance?            ==> NB
  #E. Correlation?               ==> GLMM
  #F. Non-linear patterns        ==> GAM(M) 
  #G. Wrong link function        ==> Change it 


#Standard steps to address some of these points:
#
#Fitted values
F1 <- outP$mean$mu  #<--- MCMC equivalent of fitted()
E1 <- E            #Now I can use the code from the frequentist analysis
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0)
#No outliers

#C. Zero inflation
plot(table(Owls$NCalls),
     type = "h")     
#Some zero inflation!


#D. Large variance?
dotchart(Owls$NCalls)
#That looks like a negative binomial distribution!


#E. Correlation
#Random slopes? Correlation between nests?


#F. Non-linear patterns   
#
plot(x = Owls$ArrivalTime,
     y = E1)
abline(h = 0, lty = 2)


#Two common options to deal with the overdispersion are:
# A. Add an observation level random intercept
# B. Negative binomial distribution.





##################################
#Apply a NB GLMM in JAGS.

#Use MCMC to estimate the parameters                     
#1. Bundle data
X <- model.matrix(~ FoodTreatment + ArrivalTime.std, data = Owls)                      
K <- ncol(X)  #Number of columns
head(X)

#Random effects:
Nest <- as.numeric(as.factor(Owls$Nest))
Nest
Nre <- length(unique(Owls$Nest))


#The code below is copied an pasted from the owls mixed modelling code
#Nest was changed into Hive
win.data <- list(Y     = Owls$NCalls, 
                 X     = X,
                 N     = nrow(Owls),
                 K     = K,
                 Nest  = Nest,
                 Nre   = Nre,
                 LogBS = Owls$LogBroodSize)
win.data

###################################################
# 2. JAGS modelling code

#This code was copied from the coral reef example. We only changed Site into Nest.
sink("OwlsGLMM.txt")
cat("
model{
    #1A. Priors beta and sigma
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}

    #1B. Priors random effects and sigma_Nest
    for (i in 1:Nre) { a[i] ~ dnorm(0, tau_Nest)}
    
    tau_Nest      <- 1 / (sigma_Nest * sigma_Nest)
    num           ~ dnorm(0, 0.0016)             #<----half-Cauchy(25)
    denom         ~ dnorm(0, 1)                  #<----half-Cauchy(25)
    sigma_Nest   <- abs(num / denom)             #<----half-Cauchy(25)

    #1C. Prior for size
    num.s   ~ dnorm(0, 0.0016)            #<----from squirrel exercise
    denom.s ~ dnorm(0, 1)                 #<----from squirrel exercise
    size <- abs(num.s / denom.s)          #<----from squirrel exercise


    #2. Likelihood
    for (i in 1:N) {
      Y[i] ~  dnegbin(p[i], size)
      p[i] <- size / (size + mu[i])  

      log(mu[i]) <- eta[i] 
      eta[i]     <- inprod(beta[], X[i,]) + a[Nest[i]] + LogBS[i] 
  
      #3. New stuff: Discrepancy measures 
      Exp[i] <- mu[i] 
      Var[i] <- mu[i] + mu[i] * mu[i] / size
      
      #Pearson residuals
      E[i]   <- (Y[i]  - Exp[i]) / sqrt(Var[i])    

      #Simulated NB data 
      YNew[i] ~  dnegbin(p[i], size) 

      #Pearson residuals for simulated data
      ENew[i] <- (YNew[i] - Exp[i]) / sqrt(Var[i]) 
      D[i]    <- pow(E[i], 2)                      
      DNew[i] <- pow(ENew[i], 2)   
     }          
     Fit         <- sum(D[1:N])      #Sum of squared residuals  
     FitNew      <- sum(DNew[1:N])  #Sum of squared predicted residuals
}
",fill = TRUE)
sink()
#####################################


#3. Initial values & parameters to save
inits  <- function () {
  list(
    beta  = rnorm(K, coef(M6), 0.1),
    a     = rnorm(Nre, 0, 0.1),
    num   = rnorm(1, 0, 25),    #<-- half-Cauchy(25)
    denom = rnorm(1, 0, 1),     #<-- half-Cauchy(25)
    num.s   = rnorm(1, 0, 25),  #for size
    denom.s = rnorm(1, 0, 1))   #for size
    }

params <- c("beta", "E", "a", 
            "sigma_Nest", "mu", 
            "Fit", "FitNew",
            "size")


#Don't forget to change the file name!
#4. Start JAGS
J1 <- jags(data       = win.data,
           inits      = inits,
           parameters = params,
           model      = "OwlsGLMM.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

J2  <- update(J1, n.iter = 10000, n.thin = 10)  
out <- J2$BUGSoutput



#5. Assess mixing
MyNames <- c(colnames(X), "sigma_Nest", "size")

MyBUGSChains(out, 
             c(uNames("beta", K), "sigma_Nest", "size"),
             PanelNames = MyNames)
MyBUGSACF(out, 
          c(uNames("beta", K), "sigma_Nest", "size"),
          PanelNames = MyNames)

#6. Present output
OUT1 <- MyBUGSOutput(out, 
                     c(uNames("beta", K), "sigma_Nest", "size"),
                     VarNames = MyNames)
print(OUT1, digits =5)

MyBUGSHist(out, 
           c(uNames("beta", K), "sigma_Nest", "size"),
           PanelNames = MyNames)
######################################################



######################################################
# Assess goodness of fit: Bayesion p-value:
# We have a large number of FitNew and Fit values. 
# Count how often one is bigger than the other:

mean(out$sims.list$Fit >  out$sims.list$FitNew)
#Yes...we expected that based on the frequentist results!

#Frequentist approach
E   <- out$mean$E
sum(E^2) / (nrow(Owls) - 5)


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



##################################################
#1. Apply model validation
#Extract residuals and fitted values from JAGS
E.mcmc <- out$mean$E  #Posterior mean of the residuals
F.mcmc <- out$mean$mu   #This is with random effects

#Any heterogeneity?
plot(x = F.mcmc,
     y = E.mcmc,
     xlab = "Posterior mean fitted values",
     ylab = "Posterior mean Pearson residuals")
abline(h = 0, lty = 2)     
#Ok..ish? Hmmm....


#Residual patterns?
plot(x = Owls$ArrivalTime,
     y = E.mcmc,
     xlab = "Arrival time",
     ylab = "Posterior mean Pearson residuals")
abline(h = 0, lty = 2)     

# Is there a non-linear pattern in here?
library(mgcv)
T2 <- gam(E.mcmc ~ s(ArrivalTime), data = Owls)
summary(T2)
plot(T2)
abline(h = 0, lty = 2)
# Yes...there is a non-linear pattern in here.
# We should be doing a GAMM!
# But let's go on with a GLMM for the moment


par(mfrow = c(1,1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
boxplot(E.mcmc ~ SexParent, 
        data = Owls,
        ylab = "Posterior mean residuals",
        xlab = "SexParent")
abline(h = 0, lty = 2)


par(mfrow = c(1,1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
boxplot(E.mcmc ~ FoodTreatment, 
        data = Owls,
        ylab = "Posterior mean residuals",
        xlab = "Food treatment")
abline(h = 0, lty = 2)
##############################################################



#2. Sketch the results of the model
# TASK: Write down the fitted model:
OUT1

#Now sketch the results.
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


#2. Define a grid of covariate values
range(Owls$ArrivalTime.std)
library(plyr)
MyData <- ddply(Owls, 
                .(FoodTreatment), 
                summarize,
                ArrivalTime.std = seq(min(ArrivalTime.std), 
                                      max(ArrivalTime.std),
                                      length = 25))

                       
#B. Create X matrix with expand.grid
X <- model.matrix(~ FoodTreatment + ArrivalTime.std, data = MyData)

#C. Calculate the predicted MCMC values
#   for mu
eta  <- X %*% t(beta.mcmc) + log(mean(Owls$LogBroodSize))
mu   <- exp(eta)

# Why don't we take the 2.5% and 97.5% values
# at each of the artificial covariate values?
# And plot these instead of the 3,000 lines? 
# We wrote  a small support function that does this:
	
L <- GetCIs(mu)
L

#And add all pieces in MyData
MyData2 <- cbind(MyData, L)


#Now we need some fancy plotting tools.


library(ggplot2)  
p <- ggplot()
p <- p + xlab("Arrival time") + ylab("Noises")
p <- p + theme(text = element_text(size=15))

p <- p + geom_ribbon(data = MyData2,
                     aes(x = ArrivalTime.std, 
                         ymax = up, 
                         ymin = lo),
                     fill = "red")
p <- p + geom_line(data = MyData2, 
                    aes(x = ArrivalTime.std, 
                        y = MyData2$mean, 
                        size = 1),    
                    col = "black")

p <- p + geom_point(data = Owls, 
                    aes(x = ArrivalTime.std, y = NCalls),
                    position = position_jitter(width = .02),
                    color = grey(0.7),
                    size = 2)
       
p <- p + facet_grid(. ~ FoodTreatment, 
                    scales = "fixed")
p <- p + theme(legend.position="none") 
p
########################################





########################################
#And as a 3-d graph.

library(rgl)

plot3d(x = Owls$ArrivalTime.std,
       y = Owls$LogBroodSize,
       z = Owls$NCalls,
       type = "p",
       size = 2,
       lit = FALSE,
       xlab = "ArrivalTime.std",
       ylab = "LogBroodSize",
       zlab = "NCalls")

#Now we would like to add the fitted values
#of the NB GLMM as a surface.
MyData <- ddply(Owls, 
                .(FoodTreatment), 
                summarize,
                ArrivalTime.std = seq(min(ArrivalTime.std), 
                                      max(ArrivalTime.std),
                                      length = 25))
MyData


#Calculate the fitted values
range(Owls$ArrivalTime.std)
range(Owls$LogBroodSize)
MyData <- expand.grid(ArrivalTime.std = seq(-1.59, 2.35, length = 25),
                      LogBroodSize = seq(0, 1.95, length = 25),
                      FoodTreatment = levels(Owls$FoodTreatment))
 
X    <- model.matrix(~ 1 + FoodTreatment + ArrivalTime.std, data = MyData)

eta     <- X %*% t(beta.mcmc) + MyData$LogBroodSize
mu      <- exp(eta)	        #Fitted values
L       <- GetCIs(mu)       #Take averages and CIs
MyData2 <- cbind(MyData, L) #Add to MyData 



#So..these are the two variables that we use
#to calculate our grid for:
X1 <- seq(-1.59, 2.35, length = 25)
X2 <- seq(0, 1.95, length = 25)

#And we convert the vector with expected values into
#a 25 by 25 matrix
#But...we have two planes. Extract the fitted values per FoodTreatment level
ExpY.Deprived <- MyData2$mean[MyData2$FoodTreatment == "Deprived"]
ExpY.2d.Deprived <- matrix(ExpY.Deprived, nrow = 25, ncol = 25)

ExpY.Satiated <- MyData2$mean[MyData2$FoodTreatment == "Satiated"]
ExpY.2d.Satiated <- matrix(ExpY.Satiated, nrow = 25, ncol = 25)


#Add the surface for the fitted NB model
surface3d(X1, X2, ExpY.2d.Deprived, 
          alpha = 0.6, 
          front = "lines", 
          back = "lines", 
          color = "black")

surface3d(X1, X2, ExpY.2d.Satiated, 
          alpha = 0.6, 
          front = "lines", 
          back = "lines", 
          color = "red")
          
         
#################################
#Add the 95% credible intervals!
#Deprived data
#Upper bounds
SEup.Deprived <- MyData2$up[MyData2$FoodTreatment == "Deprived"]
SEup.2d.Deprived <- matrix(SEup.Deprived, nrow = 25, ncol = 25)

#Lower bounds
SElo.Deprived <- MyData2$lo[MyData2$FoodTreatment == "Deprived"]
SElo.2d.Deprived <- matrix(SElo.Deprived, nrow = 25, ncol = 25)


#Add the surface for the fitted NB model
surface3d(X1, X2, SElo.Deprived, 
          alpha = 0.6, 
          front = "lines", 
          back = "lines", 
          color = "black")

surface3d(X1, X2, SEup.Deprived, 
          alpha = 0.6, 
          front = "lines", 
          back = "lines", 
          color = "black")


#Satiated data
#Upper bounds
SEup.Satiated <- MyData2$up[MyData2$FoodTreatment == "Satiated"]
SEup.2d.Satiated <- matrix(SEup.Satiated, nrow = 25, ncol = 25)

#Lower bounds
SElo.Satiated <- MyData2$lo[MyData2$FoodTreatment == "Satiated"]
SElo.2d.Satiated <- matrix(SElo.Satiated, nrow = 25, ncol = 25)


surface3d(X1, X2, SElo.Satiated, 
          alpha = 0.6, 
          front = "lines", 
          back = "lines", 
          color = "red")

surface3d(X1, X2, SEup.Satiated, 
          alpha = 0.6, 
          front = "lines", 
          back = "lines", 
          color = "red")



###################################################
#End of code


