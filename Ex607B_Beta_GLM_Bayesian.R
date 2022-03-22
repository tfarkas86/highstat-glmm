#    Highland Statistics Ltd.
#    www.highstat.com
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

#This exercise is based on:
# Evidence for varying social strategies across the day in chacma baboons
# Claudia Sick, Alecia J. Carter, Harry H. Marshall, Leslie A. Knapp, Torben Dabelsteen and Guy Cowlishaw
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
   

#Question 1: Subordinates should prefer to groom more dominant animals earlier in the day
#Question 2: Subordinates’ relative contributions to grooming sessions should be greater 
#            earlier in the day, especially with high-ranking partners


#What will you learn in this exercise?
# 1A. MCMC code for 2-way nested and crossed beta-GLMM



###########################################################################
#Set the working directory
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/FilesMCMC_GLMM/CDMCMCGLMGLMMCourse/Data")	
Monkeys <- read.table(file = "MonkeysV2.txt", header = TRUE)
names(Monkeys)


#########################################################################################
#Question 1: Subordinates should prefer to groom more dominant animals earlier in the day
#For this question we need to analyse only the Subordinate Grooms data:

Monkeys.sub <- Monkeys[Monkeys$SubordinateGrooms == "yes",]
MS <- Monkeys.sub  #Is shorter


#The Sick et al. (2014) paper gives a description of the experimental design, 
#but at times it is difficult to understand how exactly sampling took place.
#Here are the essential parts:

#1. Two monkey groups (large vs small). Coded as GroupSize
#2. There are multiple focal hours per day (coded as "FocalHour")
#3. In a focal hour we look for a groomer (coded as FocalGroomer)
#4. We also write down the id of the receiver monkey  (coded as Receiver)
#5. Record the time since sunrise, and the rank difference between focal groomer and receiver

#Aim: RankDifference = function(Time, Relatedness, GroupSize)

#Interactions are expected:
    #Time and relatedness  (used in paper)
    #Based on our common sense: Relatedness and groupsize  ?
########################################################################################





########################################################################################
#Understanding the dependency in the data
#Where is the correlation?
#Our first idea was:
#All observations made in the same hour
#All observations made on the same focal groomer
#All observations made on the same receiver

#This would mean:
#(1 | FocalHour) + (1 | FocalGroomer) + (1 | Receiver) 

#Sick et al use: 
#(1 | GroupSize /  FocalGroomer / FocalHour) + (1 | Receiver)
#The first one doesn't make sense (only 2 levels in GroupSize!).

#So...we better use:
#(1 | FocalGroomer / FocalHour) + (1 | Receiver)

#So..what is the difference between:
#(1 | FocalGroomer / FocalHour)
#(1 | FocalHour / FocalGroomer)
#(1 | FocalGroomer) + (1 | FocalHour)

#FocalGroomer is groomer monkey identity. 
#FocalHour is the hour in which sampling took place

#First we need to figure out how things have been coded
table(MS$FocalGroomer)
table(MS$FocalHour)
length(table(MS$FocalGroomer))
length(table(MS$FocalHour))

#59 monkeys divided over two groups
#573 hours of sampling. Note: each sampling hour is
#uniquely coded! Each moneky is uniquely coded. 
#That makes life easy!



#What does this mean: (1 | FocalGroomer / FocalHour)
#Each groomer is sampled in different hours.

#What does this mean: (1 | FocalHour / FocalGroomer)
#Per hour we sample multiple groomers.



#What do we have here?
library(lme4)
Z <- xtabs(~ FocalHour + FocalGroomer, 
           data = MS,
           drop = TRUE, sparse = TRUE)
Z
#use: head(Z)  or make your screen font smaller!

Z2 <-table(MS$FocalHour, MS$FocalGroomer)
head(Z2)
rowSums(Z2 > 0)
#So..in each hour we only follow one groomer!!

colSums(Z2 > 0)
#Most groomers are sampled in different hours.
#So we have: (1 | FocalGroomer / FocalHour)
#Make a graph of this!


#Then...what does (1 | FocalGroomer / FocalHour) do?
#Each focal groomer has a random intercept
#Each focal hour within a groomer has its own intercept
#But...focal hour is coded uniquely...hence this is the same as:
#(1 | FocalGroomer) + (1 | FocalHour)
#We will double check this in the numerical output!
############################################################################






###########################################################################
#Response variable: RankDifference  (difference between groomer and receiver)
#Covariates: Time           (time since sunrise)
#            Relatedness    (index specifying relatedness)
#            GroupSize      (small vs large)  
#
#            FocalGroomer     (monkey id who is doing the grooming)
#            FocalHour        (focal hour)
#            Receiver         (grooming partner identity)
#
#######################################################################################




#######################################################################################
#Housekeeping
library(lattice)
source(file = "/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV6.R")
library(ggplot2)
#######################################################################################





#######################################################################################
#Data exploration

#Size of the data
dim(MS)


#Outliers
MyVar <- c("RankDifference", "Time", "Relatedness")
Mydotplot(Monkeys.sub[, MyVar])
# No outliers
# Rank difference is scaled between 0 and 1!

#Collinearity
cor(MS$Time, MS$Relatedness)
#But..why would there be any correlation between
#these variables?


#Relayionships
#Create a graph with Time along the x-axis and Rank differences
#along the y-axis
p <- ggplot(data = MS, aes(x = Time , y = RankDifference))
p <- p + geom_point()  #add points
p  #Nothing yet

#Add a smoother (don't close the graph)
p <- p + geom_smooth()  #Add a smoother
p                       #Look promising!

#Split up the graph by group size
p <- p + facet_wrap(~ GroupSize)
p          # Non-linear pattern for the smaller group?






#Make a boxplot of RankDifference conditional on focal groomer 
p <- ggplot(MS, aes(y = RankDifference, x = FocalGroomer))
p <- p + geom_boxplot() 
p 
#Some groomers have consistently higher rank differences




#Is there collinearity between Time and partner ID?
p <- ggplot(MS, aes(y = Time, x = Receiver))
p <- p + geom_boxplot() 
p    #Not really
######################################################################



#Is there a relationship between Rank differnce and relatedness?
p <- ggplot(data = Monkeys.sub, aes(x = Relatedness , y = RankDifference))
p <- p + geom_point()
p <- p + geom_smooth()
p   #I don't think this is explainable
######################################################################





###################################################
#Not worked out
# print(image(xtabs(~  FocalObservation + FocalIdentity, 
            # data = MS, sparse = TRUE),
            # sub = NULL, xlab = "sample", ylab = "batch"))



#A few more graphs showing how sampling took place
dotplot(FocalGroomer ~ FocalHour, 
        data = MS,
        col = 1)


dotplot(FocalHour  ~ FocalGroomer, 
        data = MS,
        col = 1)
######################################################################







####################################################################################
#Load remaining functions and packages
source(file = "/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstatV3.R")
library(R2jags)
####################################################################################




#Step 2: Formulate JAGS modelling code
########################################################################
#What is a regression model with a beta distribution?
#See the pdf file at:  www.highstat.com/CourseMCMCGLMM/pdfs/CheetahsPage228_230.pdf
#And read Subsection 7.6.1


#The beta model is given by:
# Y ~ beta(a, b)

# E(Y)      =  a / (a + b)

#                       a * b
# var(Y)    =  ------------------------
#                 (a+b)^2 * (a + b + 1)

#Choose: 
# eta = Covariate stuff + random effects 
# logit(Pi) = eta 
# a = theta * Pi
# b = theta * (1 - Pi)

#Then we have:

# E(Y)      =  theta * Pi / (theta * Pi + theta * (1 - Pi)) = Pi    = exp(eta)  / (1 + exp(eta))

#                       a * b                      Pi * (1 - Pi)
# var(Y)    =  ------------------------    = ------------------------
#                 (a+b)^2 * (a + b + 1)             (theta + 1)


#This is similar to the binary logistic regression model...but Y can be anything between 0 and 1.
#In Bernoulli GLM Y is either 0 or 1.




#If the response variable Y is 0 < Y < 1 then the beta distribution is an option
#In the Monkey data set, we have 0 <= Y <= 1
#The response variable can take the values 0 or 1 (and for the full data set it does).
#In that case transform the variable using: (Y · (n − 1) + 0.5) / n, n is the sample size
#See: Beta regression in R. Cibrari-Neto & Zeileis. 
#     Journal of Statistical Software. April 2010, Volume 34, Issue 2.

#It can also be used if a < y < b ....see the paper.

#Rescale response variable:
N <- nrow(Monkeys.sub)
Monkeys.sub$RD.scaled <- (Monkeys.sub$RankDifference * (N - 1) + 0.5) / N
max(Monkeys.sub$RD.scaled)

#Check what the transformation has done:
plot(x = Monkeys.sub$RD.scaled, 
     y = Monkeys.sub$RankDifference,
     xlim = c(0, 1),
     ylim = c(0, 1))

Monkeys.sub$Time.std        <- Mystd(Monkeys.sub$Time)
Monkeys.sub$Relatedness.std <- Mystd(Monkeys.sub$Relatedness)


#We will first fit an ordinary beta regression model (without the random effects)
library(betareg)  #You need to install this package
B1 <- betareg(RD.scaled ~ Time.std * Relatedness.std, data = Monkeys.sub)
summary(B1)

#The model that we are fitting:
# RD.scaled_i ~ Beta(Pi_i, theta)
# E(RD.scaled_i) = Pi_i
# g(Pi_i) = Covariate stuff

# The default link is the logistic link function:
# logit(Pi_i) = Covariate stuff
# var(RD.scaled_i) = Pi_i * (Pi_i - 1) / (1 + theta)

# Coefficients (mean model with logit link):
                         # Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               0.09953    0.02382   4.179 2.92e-05 ***
# Time.std                 -0.02945    0.02386  -1.234   0.2172    
# Relatedness.std          -0.02407    0.02381  -1.011   0.3122    
# Time.std:Relatedness.std  0.05111    0.02492   2.051   0.0403 *  

# Phi coefficients (precision model with identity link):
      # Estimate Std. Error z value Pr(>|z|)    
# (phi)   4.8888     0.1842   26.54   <2e-16 ***



#Hence, the model fitted model is given by:

# RD.scaled_i ~ Beta(Pi_i, 4.8888)
# E(RD.scaled_i) = Pi_i
# logit(mu_i) = 0.099 - 0.029 * Time_i -0.02 * Relatedness_i + 0.05 * Time_i * Relatedness_i
# var(RD.scaled_i) = mu_i * (mu_i - 1) / (1 + 4.88)

#You can easily sketch this 
#########################################################################





#########################################################################
# Now let us fit the model using JAGS. First without random effects. That 
# allows us to verify that the code is correct.
# Then we will add the random effects.

#The confusing thing is that JAGS uses the following parametrisation:

# y ~ beta(a, b)
#
#           a
# E(y) = -------
#         a + b 

#                a * b
# var(y) =  -------------------
#            (a+b)^2 * (a+b+1)


#That looks different as the mean and variance that we presented for the betareg code. 
#But choose:
# a = theta * pi[i]        #(this is called shape1 in the code below)
# b = theta * (1 - pi[i])  #(this is called shape2 in the code below)

#with logit(pi) = Covariates
#And you will get exactly the same equations for the mean and the variance as presented for betareg

#To summarise...in JAGS we need a slightly different coding...but the end result is the same.
###############################################################################






###############################################################################
#JAGS coding
#Create X matrix
X <- model.matrix(~Time.std * Relatedness.std, data = Monkeys.sub)
K <- ncol(X)


#Put all required data in a list
win.data <- list(Y        = Monkeys.sub$RD.scaled,
                 N        = nrow(Monkeys.sub),
                 X        = X,
                 K        = ncol(X)
                 )

sink("BetaReg.txt")
cat("
model{
    #Priors
    for (i in 1:K) { beta[i]  ~ dnorm(0, 0.0001) }  
    theta ~ dunif(0, 10)
    
    #######################
    #Likelihood 
    for (i in 1:N){       
      #This is the 'odd' JAGS coding	
      Y[i] ~ dbeta(shape1[i], shape2[i])
      shape1[i] <- theta * pi[i]          #a
      shape2[i] <- theta * (1 - pi[i])    #b
      
      #Logistic link function
      logit(pi[i]) <- eta[i]
      eta[i]       <- inprod(beta[], X[i,]) 

      #Mean and variance, and Pearson residuals
      #ExpY[i] <- pi[i] 
      #VarY[i] <- pi[i] * (1 - pi[i])  / (theta + 1)
      #PRes[i] <- (Y[i] - ExpY[i]) / sqrt(VarY[i])

      #Discrepancy measures
      #YNew[i]     ~ dbeta(shape1[i], shape2[i])   #New data
      #PResNew[i] <- (YNew[i] - ExpY[i]) / sqrt(VarY[i])
      #D[i]       <- pow(PRes[i], 2)
      #DNew[i]    <- pow(PResNew[i], 2)
  } 
    #Fit         <- sum(D[1:N])
    #FitNew      <- sum(DNew[1:N])  
}
",fill = TRUE)
sink()





#Set the initial values for the betas and sigma
inits <- function () {
  list(
    beta  = rnorm(ncol(X), 0, 0.1),
    theta = runif(0, 10)    
    )  }

#Parameters to estimate
params <- c("beta", 
            "theta"#
            #, "PRes","Fit", "FitNew"
            )



######################################################

J0 <- jags(data = win.data,
           inits = inits,
           parameters = params,
           model.file = "BetaReg.txt",
           n.thin = 10,
           n.chains = 3,
           n.burnin = 4000,
           n.iter   = 5000)

#For fast computers:
J1<- update(J0, n.iter = 10000, n.thin = 10)  
out <- J1$BUGSoutput

#For slow computers:
out <- J0$BUGSoutput


#Assess mixing
MyBUGSChains(out, c(uNames("beta", K), "theta"))
MyBUGSACF(out, c(uNames("beta", K), "theta"))


##6. Present output
OUT1 <- MyBUGSOutput(out, c(uNames("beta", K), "theta"))
rownames(OUT1)[1:K] <- colnames(X)
print(OUT1, digits = 5)
MyBUGSHist(out, c(uNames("beta", K), "theta"))

##Is the distribution ok?
#mean(out$sims.list$FitNew > out$sims.list$Fit)


######################################################################
#Compare with betareg output (assuming betareg was installed)
B1 <- betareg(RD.scaled ~ Time.std * Relatedness.syd, data = Monkeys.sub)
summary(B1)
print(OUT1, digits = 5)
#Similar, hence our JAGS code is correct
######################################################################
######################################################################








######################################################################
######################################################################
#Add random effects to betareg model


#Covariate matrix
X <- model.matrix(~Time.std * Relatedness.std, data = Monkeys.sub)
K <- ncol(X)


#Random effects:
reFocalGroomer  <- as.numeric(as.factor(Monkeys.sub$FocalGroomer))
NumFocalGroomer <- length(unique(Monkeys.sub$FocalGroomer))

reFocalHour  <- as.numeric(as.factor(Monkeys.sub$FocalHour))
NumFocalHour <- length(unique(Monkeys.sub$FocalHour))

reReceiver  <- as.numeric(as.factor(Monkeys.sub$Receiver))
NumReceiver <- length(unique(Monkeys.sub$Receiver))


#Put all required data in a list
win.data <- list(Y            = Monkeys.sub$RD.scaled,
                 N            = nrow(Monkeys.sub),
                 X            = X,
                 K            = ncol(X),
                 reFGroomer   = reFocalGroomer,
                 NumFGroomer  = NumFocalGroomer,
                 reFHour      = reFocalHour,
                 NumFHour     = NumFocalHour,
                 reRec        = reReceiver,
                 NumReceiver  = NumReceiver
                 )
win.data

#JAGS code
sink("BetaRegmm.txt")
cat("
model{
    #Priors regression parameters, theta
    for (i in 1:K) { beta[i]  ~ dnorm(0, 0.0001) }  
    theta ~ dunif(0, 20)
    
    #Priors random effects
    for (i in 1: NumFGroomer) {a[i] ~ dnorm(0, tau.gr) }
    for (i in 1: NumFHour)    {b[i] ~ dnorm(0, tau.fh) }
    for (i in 1: NumReceiver) {c[i] ~ dnorm(0, tau.rc) }
    
    #Priors for the variances of the random effects 
    tau.gr <- 1 / (sigma.gr * sigma.gr)
    tau.fh <- 1 / (sigma.fh * sigma.fh)
    tau.rc <- 1 / (sigma.rc * sigma.rc)
    
    sigma.gr ~ dunif(0, 10)
    sigma.fh ~ dunif(0, 10)
    sigma.rc ~ dunif(0, 10)
    
    
    #######################
    #Likelihood 
    for (i in 1:N){       
      Y[i] ~ dbeta(shape1[i], shape2[i])
      shape1[i] <- theta * pi[i]        #a
      shape2[i] <- theta * (1 - pi[i])  #b
      
      logit(pi[i]) <- eta[i] + a[reFGroomer[i]] + b[reFHour[i]] + c[reRec[i]]
      eta[i]       <- inprod(beta[], X[i,]) 

      #ExpY[i] <- pi[i] 
      #VarY[i] <- pi[i] * (1 - pi[i])  / (theta + 1)
      #PRes[i] <- (Y[i] - ExpY[i]) / sqrt(VarY[i])

      ##Discrepancy measures
      #YNew[i]     ~ dbeta(shape1[i], shape2[i])   #New data
      #PResNew[i] <- (YNew[i] - ExpY[i]) / sqrt(VarY[i])
      #D[i]       <- pow(PRes[i], 2)
      #DNew[i]    <- pow(PResNew[i], 2)
  } 
    #Fit         <- sum(D[1:N])
    #FitNew      <- sum(DNew[1:N])  
}
",fill = TRUE)
sink()





#Set the initial values for the betas and sigma
inits <- function () {
  list(
    beta  = rnorm(ncol(X), 0, 0.1),
    theta = runif(0, 20),
    sigma.gr = runif(0, 10),
    sigma.fh = runif(0, 10),     
    sigma.rc = runif(0, 10)     
    )  }

#Parameters to estimate
params <- c("beta", 
            "theta", 
            #"PRes","Fit", "FitNew",
            "sigma.gr", "sigma.fh", "sigma.rc"#,
            #"ExpY"
            )



######################################################

J0 <- jags(data = win.data,
           inits = inits,
           parameters = params,
           model.file = "BetaRegmm.txt",
           n.thin = 10,
           n.chains = 3,
           n.burnin = 4000,
           n.iter   = 5000)

#Fast computer:
J1<- update(J0, n.iter = 10000, n.thin = 10)  
out <- J1$BUGSoutput

#Slow computer:
#J1<- update(J0, n.iter = 10000, n.thin = 10)  
out <- J0$BUGSoutput


#Assess mixing
MyBUGSChains(out, c(uNames("beta", K), "theta", "sigma.gr", "sigma.fh", "sigma.rc"))
MyBUGSACF(out, c(uNames("beta", K), "theta"))


##6. Present output
OUT1 <- MyBUGSOutput(out, c(uNames("beta", K), "sigma.gr", "sigma.fh", "sigma.rc", "theta"))
rownames(OUT1)[1:K] <- colnames(X)
print(OUT1, digits = 5)  

MyBUGSHist(out, c(uNames("beta", K), "theta", "sigma.gr", "sigma.fh", "sigma.rc"))


#mean(out$sims.list$FitNew > out$sims.list$Fit)
######################################################################




# # ######################################################################
# #Model validation
# E1 <- out$mean$PRes
# F1 <- out$mean$ExpY

# #Check for homogeneity
# plot(x = E1,
     # y = F1)
     
# #Check for independence
# plot(x = Monkeys.sub$Time,
     # y = E1)
# abline(h = 0)          

# plot(x = Monkeys.sub$Relatedness,
     # y = E1)
# abline(h = 0)          
# ######################################################################








######################################################################
#Sketch the fitted values

#Get the realisations of the betas
Beta.mcmc <- out$sims.list$beta 


#Create a grid of covariate values
range(Monkeys.sub$Time)
range(Monkeys.sub$Relatedness)
MyData <- expand.grid(Time        = seq(0.266, 12.22, length = 15),
                      Relatedness = seq(0, 0.75, length = 15))

#Convert the grid into a X matrix, and calculate the fitted values for each (!) MCMC iteration
X          <- model.matrix(~Time * Relatedness, data = MyData)
Betas      <- Beta.mcmc
eta        <- X %*% t(Betas)
mu         <- exp(eta) / (1 + exp(eta))
dim(mu)



#Calculate the posterior mean and 95% CI for each value on the grid
MyLinesStuff <- function(x){
   OUT <- matrix(nrow = nrow(x), ncol=4) 
	for(i in 1:nrow(x)){
	  xi <- x[i,]	
     OUT[i,3:4] <- quantile(xi, probs = c(0.025, 0.975))
     OUT[i,1] <- mean(xi)
     OUT[i,2] <- sd(xi)
	}
	colnames(OUT) <- c("mean", "se", "2.5%", "97.5%")
	OUT
}
L <- MyLinesStuff(mu)  

L    #Posterior mean, se and 95% CI for each of these covariate values

#Glue the results in L to the MyData object
MyData <- cbind(MyData,L)
MyData



#Plot the results
 #Pick on value of i
 i <- 70
 p <- wireframe(mean ~ Time + Relatedness, 
                data = MyData,
                #zlim = c(0,1),
                shade = TRUE,
                scales = list(arrows = FALSE),
                drape = TRUE, 
                colorkey = FALSE,
                screen = list(z = i, x = -60 - i /5))
 print(p)



#Add the observed data as well
mypanel <- function(x,y,z,...) {
   panel.wireframe(x,y,z,...)
   x1 <- Monkeys.sub$Time
   y1 <- Monkeys.sub$Relatedness
   z1 <- Monkeys.sub$RD.scaled
   panel.cloud(x1,y1,z1,col = 1, pch = 16,...)
}


for (i in 1:200){
 #i <- 70
 p <- wireframe(mean ~ Time + Relatedness, 
                data = MyData,
                zlim = c(0,1),
                shade = TRUE,
                scales = list(arrows = FALSE),
                drape = TRUE, 
                colorkey = FALSE,
                panel=mypanel,
                screen = list(z = i, x = -60 - i /10))
 print(p)
}





#And if you want to add the 95% credible intervals, we need to
#do some fancy coding:

#Rerun the code for the MyData again:
MyData <- expand.grid(Time        = seq(0.266, 12.22, length = 15),
                      Relatedness = seq(0, 0.75, length = 15))

#We need to copy and paste this object three times under each other:
MyData3 <- rbind(MyData, MyData, MyData)


#Now add the info for the three surfaces as an extra column:
Surfaces <- c(L[,1], L[,3], L[,4])
#The first x rows contain the posterior mean values, the
#second x rows the lower CI, and the last x rows the upper CI.
#Add this as an extra column to MyData3:
MyData3  <- cbind(MyData3, Surfaces)


#And finally make a variable for the group argument:
ID <- rep(c(1,2,3), each = nrow(MyData))
#So..this is: 1 1 1... 1 2 2 2....2 3 3 3 3.... 3
#Add this as an extra column to MyData3:
MyData3 <- cbind(MyData3, ID)
head(MyData3)

#And use the wireframe again
for (i in 1:200) {
 #i <- -40
 p <- wireframe(Surfaces ~ Time + Relatedness, 
                data = MyData3,
                group = ID,    #<--- This is the part that creates 3 surfaces
                #zlim = c(0,1),
                shade = TRUE,
                scales = list(arrows = FALSE),
                drape = TRUE, 
                colorkey = FALSE,
                screen = list(z = i, x = -60 - i /2),
                zlab = "Fit")
 print(p)
}


#You will need to do some research to see whether anyone has thought about calculating an R^2 for a beta-GLMM
#You could plot the fitted values (without the random effects) versus the observed values...I guess you get the same
#poor fit.

#END OF CODE
######################################################################
 




 