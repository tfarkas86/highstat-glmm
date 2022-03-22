#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

######################################################################

#Set the working directory and load the data
setwd("~/Dropbox/LisbonGLMM_MCMC/")
iph <- read.table(file = "IrishPh.txt",
                  header = TRUE,
                  dec = ".")
str(iph)
names(iph)
head(iph)

# Task:
# Model pH as a function of SDI
# For the moment (wrongly) ignore the other covariates.


######################################
#Load packages and support files
library(R2jags)
library(lattice)
library(ggplot2)
source("~/Dropbox/LisbonGLMM_MCMC/MCMCSupportHighstatV4.R")
######################################



#########################
#Data exploration
par(mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = iph$SDI,
     y = iph$pH,
     xlab = "SDI",
     ylab = "pH")

#########################################
#Frequentist analysis.
# This is the model that we will fit:
#  pH_i ~ N(mu_i, sigma^2)
#  E(pH_i) = mu_i
#  mu_i    = beta_1 + beta_2 * SDI_i

M1 <- lm(pH ~ SDI, 
         data = iph)
summary(M1)

#Fitted model:
#mu_i = 8.475 - 0.025 * SDI_i  
     
#Model validation
# A. Plot residuals versus fitted values
# B. Plot residuals versus each covariate in the model
# C. Plot residuals vs each covariate not in the model
# D. Check for spatial correlation
# E. Check for normality
# F. Check for influential observations

#Here we will do: A, B, C. The rest is homework.


#A. Residuals vs fitted values
E1 <- resid(M1)
F1 <- fitted(M1)

par(mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2)

#B. Residuals vs each covariate in the model
par(mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = iph$SDI, 
     y = E1,
     xlab = "SDI",
     ylab = "Residuals")
abline(h = 0, lty = 2)

#B. Residuals vs each covariate not in the model
par(mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = iph$Altitude, 
     y = E1,
     xlab = "Altitude",
     ylab = "Residuals")
abline(h = 0, lty = 2)
 
par(mar = c(5,5,2,2), cex.lab = 1.5)
boxplot(E1 ~ Forested,
        data = iph, 
        xlab = "Forested",
        ylab = "Residuals")
abline(h = 0, lty = 2)

#Home work: Check for normality, and spatial
#           independence
####################################


#Sketch fitted values
#A. Create a grid of covariate values
#B. Predict pH values for these artificial covariate values
#C. And plot the observed data, fitted values
#   and 95% CIs for the mean.

#A. Create a grid of covariate values
range(iph$SDI)
MyData <- data.frame(SDI = seq(from = 6.3, 
                               to   = 75, 
                               length = 25))

#B. Predict pH values for these artificial covariate values
P1 <- predict(M1, newdata = MyData, se = TRUE)

myfit1 <- 8.47 + 0.025 * MyData$SDI
beta <- coef(M1)
myfit2 <- beta[1] + beta[2] * MyData$SDI

X <-  model.matrix(~ 1 + SDI, data=iph) # make matrix for covariate
myfit3 <- X %*% beta # multiply coefficients by covariates to get fitted values

#Put everything in MyData (nice for ggplot2)
MyData$mu   <- P1$fit
MyData$seup <- P1$fit + 1.96 * P1$se.fit
MyData$selo <- P1$fit - 1.96 * P1$se.fit
MyData

#C. And plot the observed data, fitted values
#   and 95% CIs for the mean.
p <- ggplot()
p <- p + geom_point(data = iph, 
                    aes(y = pH, x = SDI),
                    shape = 16, 
                    size = 3)
p <- p + xlab("SDI") + ylab("pH")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_line(data = MyData, 
                   aes(x = SDI, y = mu), 
                   colour = "black")
p <- p + geom_ribbon(data = MyData, 
                     aes(x = SDI, 
                         ymax = seup, 
                         ymin = selo ),
                     alpha = 0.5)
p
##########################################





##########################################
#Now repeat the whole thing in JAGS
# We will keep the code data-set specific.


# The data exploration is the same.
# The model that we are fitting is the same.

# We will use the flowchart in Figure 2.17 
# in 'Beginner's Guide to Zero Inflated Models


# Step 1: We need to standardize the continuous covariates
#         We use the as.numeric to remove some of the text 
#         garbage produced by scale
iph$SDI.std <- scale(iph$SDI)
iph$SDI.std <- as.numeric(iph$SDI.std)


# Step 2: Bundle the data (Irish-pH specific coding)
JAGS.data <- list(pH      = iph$pH,      #Response
                  SDI.std = iph$SDI.std, #Covariate
                  N       = nrow(iph))   #Sample size
JAGS.data


###################################################
# Step 3: JAGS modelling code
sink("iphJAGS.txt")
cat(" 
model{
    #1A. Diffuse Normal priors beta and sigma
    for (i in 1:2) { beta[i] ~ dnorm(0, 0.0001)} #Precision!!

    #Diffuse uniform prior for sigma
    tau  <- 1 / (sigma * sigma)
    sigma ~ dunif(0, 20)
    
    #2. Likelihood
    for (i in 1:N) {
      pH[i]  ~ dnorm(mu[i], tau)  
      mu[i] <- beta[1] + beta[2] * SDI.std[i]   
    }
  
    #3. Calculate residuals
    for (i in 1:N) {
    	E[i] <- pH[i] - mu[i]
    }
 }           
",fill = TRUE)
sink()
#####################################


#Steps 4 and 5: Initial values & parameters to save
inits  <- function () {
  list(
    beta  = rnorm(2, 0, 0.1),
    sigma = runif(0, 20))}

#Parameters to save: beta, sigma, E and mu
params <- c("beta", "sigma", "E", "mu")


#Step 6. Run JAGS
J1 <- jags(data       = JAGS.data,
           inits      = inits,
           parameters = params,
           model      = "iphJAGS.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

J2  <- update(J1, n.iter = 10000, n.thin = 10)  
out <- J2$BUGSoutput
print(out, digits = 3)  


#Step 7: Assess mixing
MyNames <- c("Intercept", "Slope", "sigma")
MyBUGSChains(out, 
             c("beta[1]", "beta[2]", "sigma"),
             PanelNames = MyNames)


#Step 8: Look at output
MyNames <- c("Intercept", "Slope", "sigma")
OUT1 <- MyBUGSOutput(out, 
                     c("beta[1]", "beta[2]", "sigma"),
                     VarNames = MyNames)
print(OUT1, digits = 5)

#Posterior distributions
MyBUGSHist(out, 
           c("beta[1]", "beta[2]", "sigma"),
           PanelNames = MyNames)



#Step 9: Model validation
# Everything that we calculate in JAGS and
# take back comes with a whole distribution.
# So...for each observation we have a whole bunch 
# of MCMC iterations for the residuals E. The same
# for the fitted values mu.
# We can use the posterior mean values residuals and 
# posterior mean fitted values for model validation.

E  <- out$mean$E
mu <- out$mean$mu

par(mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = mu, 
     y = E,
     xlab = "Posterior mean fitted values",
     ylab = "Posterior mean residuals")
abline(h = 0, lty = 2)

#B. Residuals vs each covariate in the model
par(mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = iph$SDI, 
     y = E,
     xlab = "SDI",
     ylab = "Posterior mean residuals")
abline(h = 0, lty = 2)

#B. Residuals vs each covariate not in the model
par(mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = iph$Altitude, 
     y = E,
     xlab = "Altitude",
     ylab = "Posterior mean residuals")
abline(h = 0, lty = 2)
 
par(mar = c(5,5,2,2), cex.lab = 1.5)
boxplot(E ~ Forested,
        data = iph, 
        xlab = "Forested",
        ylab = "Posterior mean residuals")
abline(h = 0, lty = 2)

#Home work: Check for normality, and spatial
#           independence
####################################

qqnorm(E) # pretty good
hist(E, breaks=20) # yeah ok
dmat <- dist(as.matrix(iph[,2:3]), upper=TRUE, diag=TRUE)

####################################
#Step 10: Sketch the fitted values.
# Option 1: Only sketch the fitted values
# Option 2: Sketch fitted values with a 95% credible interval

# Problem for both options: SDI may not be sorted.
# So we better create some artificial SDI values on a grid
# and do predictions for these artificial values using the 
# estimated betas:

# The book website (or course website) contains a pdf that
# explain simple matrix notation. Please have a quick look 
# as we are using it here.


X <- model.matrix(~ SDI.std, data = iph)
beta <- coef(M1)
fitted(M1)  <-  X %*% beta

# A. Create artificial covariate values and convert these
#    into an X matrix.
range(iph$SDI.std)  #STANDARDIZED SCALE!!!
MyData <- data.frame(SDI.std = seq(from = -2.13, 
                                   to = 2.03, 
                                   length = 25))
X <- model.matrix(~ SDI.std, data = MyData)
X

# B. Calculate fitted values
#Under option 1:
beta <- out$mean$beta #Posterior mean betas
mu   <- X %*% beta

#Inspect results
beta
mu
dim(mu)

# C. And plot the observed data + fitted line
par(mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = iph$SDI.std, 
     y = iph$pH, 
     xlab = "Standarized SDI",
     ylab = "pH", 
     pch = 16)
lines(x = MyData$SDI.std, 
      y = mu, 
      col = 1, 
      lwd = 5)


#Under option 2:
# We have the MCMC values of the betas.
# We can calculate the fitted values for each MCMC iteration
# and use these to calculate credible intervals:

#Get the MCMC iterations for the betas
betas.mcmc <- out$sims.list$beta 
betas.mcmc
dim(betas.mcmc)   #Intercept and slope....3000 times


# Calculate the fitted values for each MCMC iteration
mu.mcmc <- X %*% t(betas.mcmc)
dim(mu.mcmc)
# This mu.mcmc contains predicted values
# for the 25 artificial SDI values....for each 
# MCMC iteration!

#We could plot these....
par(mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = iph$SDI.std, 
     y = iph$pH, 
     xlab = "Standarized SDI",
     ylab = "pH", 
     pch = 16)
for (i in 1:100){      
  lines(x = MyData$SDI.std, 
        y = mu.mcmc[,i], 
        col = i, 
        lty = 1)
}
#We only plotted fitted values for the first 100 MCMC iterations..



# Why don't we take the 2.5% and 97.5% values
# at each of the artificial covariate values?
# And plot these instead of the 3,000 lines? 
# We wrote  a small support function that does this:
	
L <- GetCIs(mu.mcmc)
L
 #GetCIs(): Support file to calculate the posterior
 #mean, se and 2.5 and 97.5 quantiles for 
 #fitted values. The input needs to be 
 #an observation - by - MCMC iteration object.

L
# Each row in L is an artificial SDI value
# The first column is the posterior mean of the MCMC iterations
# at that specific fitted value.
# The second row is the SE, the third and fourth are the 
# 2.5 and 97.5% quartiles and are used for the 95% credible intervale
# Now plot the data and draw lines for the mean,
# and 2.5% and 97.5% values:

plot(x = iph$SDI.std, 
     y = iph$pH, 
     xlab = "Standardized SDI",
     ylab = "pH", 
     pch = 16)

lines(MyData$SDI.std, L[,"mean"], col = 2, lwd = 2)
lines(MyData$SDI.std, L[,"lo"],   col = 2, lwd = 2, lty = 2)
lines(MyData$SDI.std, L[,"up"],   col = 2, lwd = 2, lty = 2)

# So..in fact you have a large number of lines lines...
# and we show the mean, 2.5% and 97.5% quartiles.
###################################################################






###################################################################
#We will now repeat the JAGS code, but make the coding more generic
#so that for another data set you only have to make minimal changes

# Step 1: We need to standardize the continuous covariates
#         We use the as.numeric to remove some of the text 
#         garbage produced by scale
iph$SDI.std <- scale(iph$SDI)
iph$SDI.std <- as.numeric(iph$SDI.std)


# Have a look at the pdf file on the website describing 
# the matrix notation

X <- model.matrix(~ SDI.std, data = iph)
head(X)
K <- ncol(X) #Number of covariates

# Step 2: Bundle the data (generic version)
JAGS.data <- list(Y = iph$pH,     #Response variable
                  X = X,          #Covariates
                  K = K,          #Number of covariates
                  N = nrow(iph))  #Sample size
JAGS.data


###################################################
# Step 3: JAGS modelling code
sink("iphJAGS2.txt")
cat(" 
model{
    #1A. Priors beta and sigma
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}

    tau  <- 1 / (sigma * sigma)
    sigma ~ dunif(0, 20)
    
    #2. Likelihood
    for (i in 1:N) {
      Y[i]  ~  dnorm(mu[i], tau)  
      mu[i] <- inprod(beta[], X[i,])   
    }
  
    #3. Calculate residuals
    for (i in 1:N) {
    	E[i] <- Y[i] - mu[i]
    }
 }           
",fill = TRUE)
sink()
#####################################


#Steps 4 and 5: Initial values & parameters to save
inits  <- function () {
  list(
    beta  = rnorm(K, 0, 0.1),
    sigma = runif(0, 20))}

#Parameters to save: beta, sigma, E and mu
params <- c("beta", "sigma", "E", "mu")


#Step 6. Run JAGS
JG1 <- jags(data       = JAGS.data,
           inits      = inits,
           parameters = params,
           model      = "iphJAGS2.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

JG2  <- update(JG1, n.iter = 10000, n.thin = 10)  
outG <- JG2$BUGSoutput
print(outG, digits = 3)  


#Step 7: Assess mixing
MyNames <- c("Intercept", "Slope", "sigma")
uNames("beta",K) #Simple function to get a vector with names of the
                 #form:  beta[1], beta[2], ..., beta[K] 
MyBUGSChains(out, 
             c(uNames("beta", K), "sigma"),
             PanelNames = MyNames)


#Step 8: Look at output
MyNames <- c("Intercept", "Slope", "sigma")
OUT1 <- MyBUGSOutput(out, 
                     c(uNames("beta", K), "sigma"),
                     VarNames = MyNames)
print(OUT1, digits = 5)

#Posterior distributions
MyBUGSHist(out, 
           c(uNames("beta", K), "sigma"),
           PanelNames = MyNames)
##############################################






##############################################
# If you want to make these graphs yourself, use:
out$sims.matrix
colnames(out$sims.matrix)
out$sims.matrix[,"beta[1]"]
out$sims.matrix[,c("beta[1]", "beta[2]")]

out$sims.array
out$sims.array[,,"beta[1]"]
out$sims.array[,,"beta[2]"]

out$sims.list
names(out$sims.list)
out$sims.list$beta
