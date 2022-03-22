#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
# Squirrel example
# Set the working directory and import the data

setwd("~/Dropbox/LisbonGLMM_MCMC/")
SQ <- read.table(file   = "RedSquirrelsV2.txt", 
                 header = TRUE, 
                 dec    = ".")

str(SQ)
names(SQ)
###################################################################




###################################################################
# Load packages and library files
library(lattice)  #Needed for multi-panel graphs
source(file = "HighstatLibV9.R")  
# This file can be downloaded from the course website. 
# Save it in the same directory as your working directory.
# That is easier.
##################################################################


##################################################################
# Data exploration
# Outliers
MyVar <- c("SqCones", "Ntrees", "DBH", 
           "TreeHeight", "CanopyCover")

Mydotplot(SQ[,MyVar])
# Ouch....1 large DBH value!
# A few small Canopy values!
# A few large Ntrees values.
#   There may be a patch with large DBH and NTrees value!


#Remove the value with extreme large DBH value
SQ2 <- subset(SQ, DBH < 0.6)
dim(SQ)
dim(SQ2)  #That's a small data set!
          #Remember the rule: 10-15-20 observations per 
          #                   parameter (= beta)


#Collinearity
MyVar <- c("Ntrees", "DBH", 
           "TreeHeight", 
           "CanopyCover")

pairs(SQ2[, MyVar], 
      lower.panel = panel.cor)

corvif(SQ2[,MyVar])
# If you are not familiar with VIF values, see our 2010 paper.
# In principle that is ok!!
# But the outliers may be doing funny things!


# Relationships
MyVar <- c("Ntrees", "DBH", 
           "TreeHeight", 
           "CanopyCover")
Myxyplot(SQ2, MyVar, "SqCones", 
         MyYlab = "Number of stripped cones")
# Not much going on!


# Zero inflation?
plot(table(SQ2$SqCones))
sum(SQ2$SqCones == 0)  #Number of zeros
100 * sum(SQ2$SqCones == 0) / nrow(SQ2)  #% of zeros

# Based on our experience, you should raise your eyebrow
# if you have nore than 20% of zeros. But this may be
# different for your data set.


# Conclusion data exploration: 
#  -1 outlier in DBH
#  -Small amount of collinearity
#  -Weak patterns between Y and X
#  -Collinearity?
#  -Potentially NB distribution for SqCones
#####################################################





#####################################################
#Start frequentist analysis

######################################
# Step 1: Fit a Poisson GLM:

# SqCones_i ~ Poisson(mu_i)
# E(SqCones_i)   = mu_i
# var(SqCones_i) = mu_i
# log(mu_i) = alpha + beta_1 * Ntrees.std_i + .... 
#                     beta_3 * CanopyCover_std_i

M1 <- glm(SqCones ~ Ntrees +  TreeHeight + CanopyCover,
          family = "poisson",
          data = SQ2)


# Model validation
# Assess overdispersion using Pearson residuals
# What are Pearson residuals?

#   Y - E(Y)         Y - mu
# -----------  =  ------------ 
# sqrt(var(Y))      sqrt(mu)

# In R code:
E1 <- resid(M1, type = "pearson")
N  <- nrow(SQ2)
p  <- length(coef(M1))
sum(E1^2) / (N - p)   # Should be 1-ish. At the end of this
                      # solution file we will do a small simulation
                      # study to get a feel for the amount of variation
                      # that you may expect.
# Here it is: 13.68
# That is way too large!


# Why do we have overdispersion?
# A. Outliers Y?           ==> Remove them?
# B. Missing covariates?   ==> Add them (or add a latent variable)
# C. Missing interactions? ==> Add them
# D. Zero inflation?       ==> ZIP
# E. Large variance?       ==> NB or Generalized Poisson (Ntzoufras 2009)
# F. Correlation?          ==> GLMM
# G. Non-linear patterns   ==> GAM 
# H. Wrong link function   ==> Change it 

# Your task is to find the culprit. If you pick the wrong one,
# then you may end up with biased parameters.
# How do you figure out which one to pick?
#  -Know your data
#  -Model validation
#  -Data exploration


# Standard steps to address some of these points:
F1 <- fitted(M1)

par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x = F1, 
     y = E1, 
     xlab = "Fitted values",
     ylab = "Pearson residuals", 
     cex.lab = 1.5)
abline(h = 0, lty = 2)
     
#A. Outliers
plot(cooks.distance(M1), 
     type = "h",
     xlab = "Observation", 
     ylab = "Cook distance",
     cex.lab =  1.5)

# Patterns
plot(y = SQ2$SqCones, 
     x = F1,
     xlab = "Fitted values",
     ylab = "Observed data",
     cex.lab = 1.5,
     xlim = c(0,60), 
     ylim = c(0,60) )
abline(coef = c(0, 1), lty = 2)   

#D. Zero inflation
100 * sum(SQ2$SqCone==0) / nrow(SQ2)  #=10%


#G. Non-linear patterns
MyVar <- c("Ntrees", "DBH", 
           "TreeHeight", 
           "CanopyCover")

SQ2$E1 <- E1    
Myxyplot(SQ2, MyVar, "E1", 
         MyXlab = "Covariate",
         MyYlab = "Pearson residuals")
T1 <- gam(E1 ~ s(DBH), data = SQ2) #Or do a GAM with each covariate
summary(T1)                        #If this GAM shows that there is a 
                                   #significant effect of DBH, then 
                                   #the model is not good and needs to
                                   #be improved.

T2 <- gam(E1 ~ s(Ntrees), data = SQ2) #Or do a GAM with each covariate
summary(T2)

T3 <- gam(E1 ~ s(TreeHeight), data = SQ2) #Or do a GAM with each covariate
summary(T3)

T4 <- gam(E1 ~ s(CanopyCover), data = SQ2) #Or do a GAM with each covariate
summary(T4)

#Etc.

     
# We cannot pinpoint any of the common causes for
# overdispersion. The variation in the data is relatively
# large...try a negative binomial?     
###############################################



###############################################
# Negative binomial GLM
# This is the model that we will fit:

# SqCones_i      ~ NB(mu_i, theta)
# E(SqCones_i)   = mu_i
# var(SqCones_i) = mu_i + mu_i^2 / theta
# log(mu_i)      = Intercept + beta_1 * NTrees_i + ...... 


# Don't mix these up in your head:
#  -dispersion statistic: sum(E^2) / (N-p)
#  -overdispersion/underdispersion: compare dispersion 
#                                   statistic > 1  or  < 1
#  -theta: variance parameter (sister of sigma^2 from lm)


library(MASS)
M2 <- glm.nb(SqCones ~ Ntrees + TreeHeight + CanopyCover,
             data = SQ2)

# Dispersion statistic:
E2 <- resid(M2, type = "pearson")
N  <- nrow(SQ2)
p  <- length(coef(M2)) + 1  # '+1' is due to theta
sum(E2^2) / (N - p)
# Minor underdispersion..so we can inspect the numerical output


#Numerical output of the NB GLM:
summary(M2)


#What is the fitted model?
# SqCones_i      ~ NB(mu_i, 1.079)
# E(SqCones_i)   = mu_i
# var(SqCones_i) = mu_i + mu_i^2 / 1.079

#              -4.18 + 0.02 * NTrees_i + ...+ 0.06 * CanopyCover_i 
# mu_i      = e 



#Generalized R^2:

#                   Null deviance - residual deviance
#Generalized R^2 = ------------------------------------
#                           Null deviance
(77.262 - 59.073) / 77.262
#Small-ish


#Step 2: Is everything significant?
drop1(M2, test = "Chi")

#What is drop one doing?
#Answer: 1. First it fits the model with all 3 covariates
#        2. Then it fits 3 sub-models. In each sub-model 
#           one covariate is dropped. 
#        3. It then compares the deviances of the full model
#           and a submodel. The difference in deviances is
#           Chi-square distributed. The df is the difference in
#           number of parameters. 

# # drop1compares models holding theta constant for each model
# # If instead you compare seperately fit models, they will have different theta
# # and results will change



#Some of the terms are not significant
#Option 1: Leave the model as it is.
#Option 2: Drop covariates using AIC or based on p-values (considered bad).
step(M2)  #Backwards selection using AIC




# Based on the AIC (+ difference < 2 go for more simple model),
# this seems to be the optimal model:
M3 <- glm.nb(SqCones ~ CanopyCover,
             data = SQ2)
summary(M3)



######################################
# Model validation of the optimal NB GLM


# Plot residuals vs fitted values
F3 <- fitted(M3)
E3 <- resid(M3, type = "pearson")
plot(x = F3, 
     y = E3,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)
     
     
# Plot the residuals vs each covariate     
SQ2$E3 <- E3
Myxyplot(SQ2, MyVar, "E3")
# Don't get too fusy about points at the edges.
# We can't really fault the model. The small sample 
# size does not make it easy neither.
###########################################################





###########################################################
#Model interpretation of the NB GLM
#Sketch the fitted values

par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = SQ2$CanopyCover, 
     y = SQ2$SqCones,
     xlab = "Canopy cover",
     ylab = "Number of stripped cones",
     cex.lab = 2,
     pch = 16, 
     type = "p")

# Create an artifical grid of covariate values
MyData <- data.frame(CanopyCover = seq(from  = min(SQ2$CanopyCover),
                                       to    = max(SQ2$CanopyCover), 
                                       length = 25))
                                       
# Predict the expected squirrel values                                       
P1 <- predict(M3, newdata = MyData, type = "link", se = TRUE)
lines(x = MyData$CanopyCover, 
      y = exp(P1$fit), 
      lwd = 3)
lines(x = MyData$CanopyCover, 
      y = exp(P1$fit + 1.96 * P1$se.fit), 
      lwd = 3, 
      lty = 2)
lines(x = MyData$CanopyCover, 
      y = exp(P1$fit - 1.96 * P1$se.fit), 
      lwd = 3, 
      lty = 2)
     
#############################################






############################################
# Let's visualise what exactly the difference between
# a Poisson GLM and a NB GLM is.

M4.pois <- glm(SqCones ~ CanopyCover,
               family = poisson,
               data = SQ2)

M4.poisq <- glm(SqCones ~ CanopyCover,
               family = quasipoisson,
               data = SQ2)

M4.nb <- glm.nb(SqCones ~ CanopyCover,
                data = SQ2)

summary(M4.pois)
summary(M4.poisq)
summary(M4.nb)  
  
  
##################################  
# Calculate fitted values.
# The code below assumes that you are familiar
# with some basic matrix notation. See the online pdf file.
# Fitted values aof the NB GLM are given by: mu = exp(X * beta)

X    <- model.matrix(M4.nb)
beta <- coef(M4.nb)
MyFittedValues <- exp(X %*% beta)

#Test it with the fitted() function:
MyFittedValues - fitted(M4.nb)  #Should be 0
##################################




# Create an artificial covariate
MyData <- data.frame(CanopyCover = seq(from = min(SQ2$CanopyCover),
                                       to   = max(SQ2$CanopyCover),
                                       length = 25))

#Create a design matrix X for **THESE** 25 values
X <- model.matrix(~ CanopyCover, 
                  data = MyData)

#Calcuate the fitted Poisson and NB values
MyData$mu.pois <- exp(X %*% coef(M4.pois))
MyData$mu.nb   <- exp(X %*% coef(M4.nb))



# We are now going to make a 3-d scatterplot.
# The code is not so nice on the eyes. Don't try
# to fully understand it.


library(scatterplot3d) 
library(VGAM)

# Ensure these ones are installed

# VGAM overwrites a function from mgcv.
# So you better close and restart R AFTER this exercise



# COPY THE R CODE FROM HERE .............

# Three variable for 3-d scatterplot
x      <- MyData$CanopyCover
y.pois <- MyData$mu.pois
y.nb   <- MyData$mu.nb
z      <- 0 * x


par(mfrow = c(1, 1))
rr<- scatterplot3d(SQ2$CanopyCover, 
                   SQ2$SqCones,
                   rep(0, nrow(SQ2)), 
                   highlight.3d = FALSE, 
                   col.axis = "black",
                   col.grid = "black", 
                   pch = 20,
                   zlim = c(0, 1),
                   ylim = c(0, max(SQ2$SqCones)),
                   type="p",
                   grid = FALSE,
                   box = TRUE,
                   cex.lab = 1.5,
                   xlab = "CanopyCover",
                   ylab = "Possible values",
                   zlab = "Probability")

#Poisson fitted values in red
rr$points3d(x, y.pois, z, 
            lwd = 3, 
            col = 2,
            type = "l")

#Negative binomial fitted values in green
rr$points3d(x, y.nb, z, 
            lwd = 3, 
            col = 3,
            type = "l")

#Pick 5 values along the x axis to plot the density curves.
MyDatai <- data.frame(CanopyCover = seq(from = quantile(MyData$CanopyCover, 0.15),
                                        to   = 0.99 * max(MyData$CanopyCover),
                                        length = 10))


K <- length(MyDatai$CanopyCover)
for (i in 1:K){
  Xi   <- model.matrix(~ CanopyCover, data = MyDatai)[i,]
  yseq <- round(seq(0, 60, by = 1))

  #Density curve Poisson GLM M3
  muPois.i <- exp(Xi %*% coef(M4.pois))
  zi       <- dpois(yseq, lambda = muPois.i)

  #Density curve NB
  muNB.i <- exp(Xi %*% coef(M4.nb))
  zi.nb <- dnbinom(yseq, size = M4.nb$theta, mu = muNB.i)
   
  #Dotted line 
  rb = cbind(MyDatai$CanopyCover[i], yseq, 0)
  rr$points3d(rb, col = 1, type = "l", pch = ".", lty = 2)
 
  #Density curves Poisson
  rb = cbind(MyDatai$CanopyCover[i]- 0.05, yseq, zi)
  rr$points3d(rb, col = 2, type = "h", pch = ".", lwd = 1)
  
  #Density curves NB
  rb = cbind(MyDatai$CanopyCover[i], yseq, zi.nb)
  rr$points3d(rb, col = 3, type = "h", pch = ".", lwd = 2)
  
 }
 
# ALL THE WAY UP TO HERE..AND RUN IN R.

# Conclusions:
#  - The Poisson and NB fitted values are nearly the same.
#  - The difference is in the range of allowable values.
#     -Red is Poisson
#     -Green is NB
############################################ 




 
 
 
 
############################################
#Same graph using another package
library(rgl)

#First we plot the observed data
plot3d(x = SQ2$CanopyCover,
       y = SQ2$SqCones,
       z = rep(0, nrow(SQ2)),
       type = "s",
       size = 1,
       lit = FALSE,
       zlim = c(0, 0.5),
       xlab = "CanopyCover",
       ylab = "Stripped cones",
       zlab = "Probability")

#You can spin this graph!
#Don't close it!

#Add the Poisson fitted values
lines3d(x, y.pois, z, 
        col = 2,
        lwd = 3)

#Add the NB fitted values
lines3d(x, y.nb, z, 
        col = 2,
        lwd = 3)

#And add the density curves. Copy and paste the whole block of
#code from here......
interleave <- function(v1, v2) as.vector(rbind(v1, v2))

for (i in 1:K){
  Xi   <- model.matrix(~ CanopyCover, data = MyDatai)[i,]
  yseq <- round(seq(0, 60, by = 1))

  #Density curve Poisson GLM M3
  muPois.i <- exp(Xi %*% coef(M4.pois))
  zi       <- dpois(yseq, lambda = muPois.i)

  
  #Density curve NB
  muNB.i <- exp(Xi %*% coef(M4.nb))
  zi.nb  <- dnbinom(yseq, size = M4.nb$theta, mu = muNB.i)
   
  #Dotted line 
  rb = cbind(MyDatai$CanopyCover[i], yseq, 0)
  lines3d(rb[,1], rb[,2], rb[,3], 
          col = 1,
          lty = 2)

  #Density curves Poisson
  rb = cbind(MyDatai$CanopyCover[i]- 0.05, yseq, zi)
  #lines3d(rb[,1], rb[,2], rb[,3], 
  #        col = 2,
  #        lwd = 1)
  segments3d(interleave(rb[,1], rb[,1]),
             interleave(rb[,2], rb[,2]),
             interleave(rep(0, nrow(rb)), rb[,3]),
             alpha = 0.9,
             col = "green")
  
  #Density curves NB
  rb = cbind(MyDatai$CanopyCover[i], yseq, zi.nb)
  #lines3d(rb[,1], rb[,2], rb[,3], 
  #        col = 3,
  #        lwd = 1)
 
 
  segments3d(interleave(rb[,1], rb[,1]),
             interleave(rb[,2], rb[,2]),
             interleave(rep(0, nrow(rb)), rb[,3]),
             alpha = 0.9,
             col = "purple")
 }
 
# All the way up to here!!
############################################








############################################
# Simulate from the model
# This is often ignored....but can be very important.
# You should do this as part of the model validation,
# or model explanation.


# The R code below is getting on the ugly side. There
# is no need to try and fully understand the R coding.
# Just run the code and try to grasp the idea and
# motivation behind the graphs.



# Let's simulate from the (overdispersed) Poisson
# GLM first.

# This was the full Poisson GLM:
M1 <- glm(SqCones ~ Ntrees +  TreeHeight + CanopyCover,
          family = "poisson",
          data = SQ2)

# So the fitted values are:
mu <- fitted(M1)

# We can also calculate these manually:
X <- model.matrix(~Ntrees +  TreeHeight + CanopyCover,
                   data = SQ2)
beta <- coef(M1)                   
mu   <- exp(X %*% beta)

#Check our results:
fitted(M1) - mu    #Bingo



# Now we do a simulation study.
# We simulate Poisson data by simulating
# data from the rpois function where the mean 
# is given by the fitted values of the model. 
N    <- nrow(SQ2)                      # Sample size
NSim <- 1000                           # Number of simulated data sets
YPois <- matrix(nrow = N, ncol = NSim) # Create space for the 1,000 data sets
for (i in 1:NSim) {
	YPois[,i] <- rpois(N, lambda = mu)
}
# Now we have 1000 simulated data sets from the model.
# What shall we do with these simulated data sets?




# We could calculate the number of zeros in each of the 1,000
# data sets.
zeros <- vector(length = 1000)
for(i in 1:NSim){
	zeros[i] <- sum(YPois[,i] == 0)
}

table(zeros)
#  0   1   2   3 
#536 386  72   6 
# Your numbers may differ due to random sampling.

# From the 1,000 simulated data sets, in 536 simulated
# data sets we had 0 zeros. In 386 simulated data sets
# we had 1 zero, in 72 data sets we had 2 zeros, and in 6
# data sets we had 3 zeros.
# Just type 
YPois[,1]
YPois[,2]
YPois[,3]
#etc



#Let's plot this as a table
plot(table(zeros), 
     #axes = FALSE,
     xlab = "How often do we have 0, 1, 2, 3, etc. number of zeros",
     ylab = "Number of zeros in 1000 simulated data sets",
     xlim = c(0, 5),
     main = "Simulation results")
points(x = sum(SQ2$SqCones == 0), 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#The red dot is the number of zeros in the original data set.
#
# TASK: What does this tell you?
# Answer: The data simulated from the Poisson model
#         does not contain enough zeros.
#########################################






#########################################
# What else can we do with the 1,000 simulated data sets?
# We can also calculate the dispersion statistic for each
# simulated data set

#RUN FROM HERE........
p <- length(coef(M1))    #Number of regression parameters
Dispersion <- vector(length = NSim) #Create space

# Calculate the dispersion for each simulated data set
for(i in 1:NSim){
	e1 <- (YPois[,i] - mu) /sqrt(mu)
	Dispersion[i] <- sum(e1^2) / (N - p)
}

# Plot this as a table
hist(Dispersion,
     xlab = "Dispersion",
     ylab = "Frequency",
     xlim = c(0, 15))

# And visualize the dispersion for the original Poisson GLM
DispersionModel <- sum(E1^2) / (N-p)
points(x = DispersionModel, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
# The red dot is the number of zeros in the original data set.
# We have serious overdispersion!

#UP TO HERE....
#########################################################################





#########################################################################
# Let's simulate from the NB GLM.
# This was the model
M5 <- glm.nb(SqCones ~ CanopyCover,
             data = SQ2)
summary(M5)

# Get the Pearson residuals
E5 <- resid(M5, type = "pearson")

# The fitted values and theta (or k) are:
mu    <- fitted(M5)
theta <- summary(M5)$theta



# Now we do a simulation study:
N    <- nrow(SQ2)   #Sample size
NSim <- 1000        #Number of simulated data sets
YNB <- matrix(nrow = N, ncol = NSim) #Create space for the 1,000 data sets
for (i in 1:NSim) {
	YNB[,i] <- rnegbin(N, mu = mu, theta = theta)
}
# Now we have 1000 simulated data sets from the model.


# What shall we do with these simulated data sets?
# We could calculate the number of zeros in each of the 1,000
# data sets.
zeros <- vector(length = 1000)
for(i in 1:NSim){
	zeros[i] <- sum(YNB[,i] == 0)
}

table(zeros)
# From the 1,000 simulated data sets, in 5 simulated
# data sets we had 0 zeros. In 67 simulated data sets
# we had 1 zero, in 124 data sets we had 2 zeros, etc.



#Let's plot this as a table
plot(table(zeros), 
     #axes = FALSE,
     xlab = "Percentage of zeros",
     ylab = "Frequency",
     xlim = c(0, 10))
points(x = sum(SQ2$SqCones == 0), 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
# The red dot is the number of zeros in the original data set.

# Question: What does this tell you?
# Answer:   The data simulated from the NB model
#           does contain enough zeros.



# What else can we do with the 1,000 simulated data sets?
# We can also calculate the dispersion statistic for each
# simulated data set
p <- length(coef(M5)) + 1  #'+1' is for the theta parameter
Dispersion <- vector(length = NSim)
for(i in 1:NSim){
	ei <- (YNB[,i] - mu) / sqrt(mu + mu^2 / theta)
	Dispersion[i] <- sum(ei^2) / (N - p)
}

#Let's plot this as a table
hist(Dispersion,
     xlab = "Dispersion",
     ylab = "Frequency",
     xlim = c(0, 3.5), breaks=40)

DispersionModel <- sum(E5^2) / (N-p)
points(x = DispersionModel , y = 0, pch = 16, cex = 5, col = 2)
abline(v=DispersionModel, lwd=5, col=2)

# The red dot is the number of zeros in the original data set.
# That is quite a large variation in dispersion statistics.
# This can be due to the very small sample size.



