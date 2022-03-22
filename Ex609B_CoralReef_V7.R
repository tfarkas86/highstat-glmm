
#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
#Coral reef exercise
#These data are taken from:
#Caribbean-wide decline in carbonate production threatens 
#coral reef growth#Perry et al. (2012)
#NATURE COMMUNICATIONS | 4:1402 | DOI: 10.1038/ncomms2409 | www.nature.com/naturecommunications


#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Books/BGS/GAMM/Data/ReefData")
CO <- read.table("CoralData.txt", header = TRUE)

names(CO)
str(CO)
# 'data.frame':	101 obs. of  10 variables:
 # $ Country         : Factor w/ 4 levels "Bahamas","Belize",..: 1 1  ...
 # $ Site            : Factor w/ 14 levels "BAR","CAI","CAL",..: 5 5  ...
 # $ Transect        : int  1 2 3 1 2 3 1 2 3 1 ...
 # $ Depth           : int  19 19 19 17 17 17 20 20 20 5 ...
 # $ Habitat         : Factor w/ 5 levels "APZ","FRMZ","FRS",..: 4 4  ...
 # $ LCC             : num  7.88 9.61 5.69 5.11 4.09 ...
 # $ Gross_production: num  2.04 1.94 1.25 1.1 1.06 1.33 1.8 1.37 2. ...
 # $ Gross_erosion   : num  1.49 1.51 1.54 1.7 1.71 1.64 1.2 1.17  ...
 # $ Net_production  : num  0.56 0.43 -0.29 -0.6 -0.65 -0.31 0.6  ...
 # $ Accretion_rate  : num  0.57 0.49 0.05 -0.21 -0.23 -0.01 0.57  ... 


###################################################################
#Load packages and library files
library(lattice)  #Needed for multi-panel graphs
library(R2jags)
library(nlme)
library(mgcv)

source(file = "/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV7.R")  
source(file = "/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstatV2.R")
##################################################################

#The variables:
# "Country"          "Site"             "Transect"         "Depth"
# "Habitat"          "LCC"              "Gross_production" "Gross_erosion"
# "Net_production"   "Accretion_rate"

# Aim of the study:
# Perry et al. (2013) sampled reef carbonate 
# production (CaCO3 m–2 yr–1) at 101 transects 
# in 19 coral reefs occurring in four Caribbean 
# countries (Bahamas, Belize, Bonaire, and Grand Cayman). 
# Surveys were conducted within a range of common 
# Caribbean reef habitats: nearshore hardgrounds, 
# Acropora palmata habitats, Montastraea spur-and-groove zones, 
# fore-reef slopes, and deeper (18–20 m) shelf-edge Montastraea reefs. 
# See Perry et al. (2013) for further details.


# The response variable is Net_production, 
# and we have one observation within a transect. 
# Hence, transect is the statistical sampling unit. 
# Potential explanatory variables are Country, 
# Site (a group of transects), Depth, Habitat, and 
# LCC (life coral cover %). 

# Task:
# Model net carbonate production as a function of:
# LCC ( Live coral cover)
# Habitat
# Random effect site nested within a Country
###################################################




###################################################
#Housekeeping
CO$fCountry <- factor(CO$Country)
CO$fSite    <- factor(CO$Site)
CO$fHabitat <- factor(CO$Habitat)

#I'm lazy....this is the response variable
CO$G        <- CO$Net_production
###################################################



###################################################
#Data exploration
#Outliers
MyVar <- c("G", "LCC")
Mydotplot(CO[,MyVar])


table(CO$fCountry)
table(CO$fSite)
table(CO$fHabitat)

#4 countries
#14 uniquely defined sites


#Relationships
boxplot(G ~ fSite, data = CO, xlab = "Site")
boxplot(G ~ fCountry, data = CO, xlab = "Country")
boxplot(G ~ fHabitat, data = CO, xlab = "Habitat")

MyYlab <- expression(paste("Net carbonate production ", 
                "(kg CaCO"[3],
                " m"^"-2",
                "year" ^"-1",")"))
par(mar = c(5,6,2,2), cex.lab = 1.5)
plot(x = CO$LCC, 
     y = CO$G,
     pch = 16,
     cex = 1,
     xlab = "LCC (%)",
     ylab = MyYlab)
abline(h = 0, lty = 2)


#############
#Collinearity
boxplot(LCC ~ fSite, data = CO, xlab = "Site", ylab = "LCC (%)")
boxplot(LCC ~ fCountry, data = CO, xlab = "Country", ylab = "LCC (%)")

par(cex.lab = 1.5, mar = c(5,5,2,2))
boxplot(LCC ~ fHabitat, data = CO, xlab = "Habitat", ylab = "LCC (%)")

#Some strong collinearity going on here!

#There is more trouble:
table(CO$fCountry, CO$fHabitat)
#Use them both an you may get trouble!
#Some levels of country cannot be discriminated from some habitat levels.
#Options: 
# 1 Remove Bahamas and SEMR 
# 2. Drop one covariate and be careful with the interpretation of the results


# Some data phishing indicates that country is not significant...but 
# habitat level is!
# Try presenting that in a paper!


#Data exploration indicates that there are no obvious 
#outliers; there is a positive, potentially nonlinear 
#LCC effect; there is heterogeneity; and there is a 
#certain degree of collinearity between LCC and the 
#categorical covariates. There is also trouble with country and habitat.
##########################################################





##########################################################
#Start of Bayesian analysis. 

#Obviously...we should start with G...apply JAGS...discover
#heterogeneity...and then decide to go for a transformation... 
#But I am saving you some time here:
CO$LogG <- log(CO$G + 3)  #I really prefer a GAMM here! 


#With start with the same model as the frequency analysis
CO$LCC.std <- (CO$LCC - mean(CO$LCC)) / sd(CO$LCC)
#A common mistake we see is that people make mistakes with the brackets!!!
 
                   
                       
#Step 1: Fit the model in JAGS.

                  
#Reference model for MCMC 
M0 <- lme(LogG ~ LCC.std + fHabitat,
                       random =~ 1 |fSite,
                       data = CO, method = "REML")
#By design of the study we use site as random effect.
                
                  
##################################################                     
#Use MCMC to estimate the parameters                     
#1. Bundle data
X <- model.matrix(~ LCC.std + fHabitat, data = CO)
#X <- model.matrix(~ LCC.std, data = CO)
#X <- model.matrix(~ fHabitat, data = CO)
    
    
    
                   
K <- ncol(X)  #Number of columns
head(X)

#Random effects:
Site <- as.numeric(as.factor(CO$fSite))
Site
Nre <- length(unique(CO$fSite))

win.data <- list(Y    = CO$LogG,  #Response
                 X    = X,        #Covariates
                 N    = nrow(CO), #Sample size
                 K    = K,        #Number of betas
                 Site = Site,     #random effect index
                 Nre  = Nre)      #number of random effects
win.data

###################################################
# 2. JAGS modelling code

sink("coralmixedmodel.txt")
cat("
model{
    #1A. Priors beta and sigma
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}
    
    tau <- 1 / (sigma * sigma)
    sigma ~ dunif(0, 10)

    #1B. Priors random effects and sigma_Site
    for (i in 1:Nre) { a[i] ~ dnorm(0, tau_Site)}
    
    tau_Site <- 1 / (sigma_Site * sigma_Site)
    sigma_Site ~ dunif(0, 10)

    #2. Likelihood
    for (i in 1:N) {
      Y[i]   ~ dnorm(mu[i], tau)   
      mu[i]  <- eta[i] + a[Site[i]]
      eta[i] <- inprod(beta[], X[i,]) 
 
      #3. Discrepancy measures 
      Exp[i] <- mu[i]      #Expected value of Y
      Var[i] <- sigma^2    #Variance of Y 
      E[i]   <- (Y[i]  - Exp[i]) / sqrt(Var[i])    
      #Normalized residuals
             
     }          
}
",fill = TRUE)
sink()
#####################################
#

#3. Initial values & parameters to save
inits  <- function () {
  list(
    beta       = rnorm(K, 0, 0.1),
    sigma      = runif(1, 0, 10),
    a          = rnorm(Nre, 0, 0.1),
    sigma_Site = runif(1, 0, 10))  }


params <- c("beta", "sigma", "E", "a", 
            "sigma_Site", "mu", "eta")


#4. Start JAGS
J1 <- jags(data       = win.data,
           inits      = inits,
           parameters = params,
           model      = "coralmixedmodel.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

J2 <- update(J1, n.iter = 10000, n.thin = 10)  
print(J2, digits = 3)

#5. Assess mixing
out <- J2$BUGSoutput

MyNames <- c(colnames(X), "Sigma", "Sigma site")

MyBUGSChains(out, 
             c(uNames("beta", K), "sigma", "sigma_Site"),
             PanelNames = MyNames)
             
MyBUGSACF(out, 
          c(uNames("beta", K), "sigma", "sigma_Site"),
          PanelNames = MyNames)
          

#6. Present output
OUT1 <- MyBUGSOutput(out, 
                     c(uNames("beta", K), "sigma", "sigma_Site"),
                     VarNames = MyNames)
rownames(OUT1)[1:K] <- colnames(X)
print(OUT1, digits =5)

MyBUGSHist(out, 
           c(uNames("beta", K), "sigma", "sigma_Site"),
           PanelNames = MyNames)


OUT1
summary(M0)  #That is all good!
######################################################




#Model validation
E   <- out$mean$E
Fit <- out$mean$mu 

par(mfrow = c(2,2), mar = c(5,5,2,2))    
plot(x = Fit, 
     y = E, 
     xlab = "Fitted values", 
     ylab = "Residuals")
abline(h = 0, lty = 2)

plot(x = CO$LCC, 
     y = E, 
     xlab = "LCC", 
     ylab = "Residuals")
abline(h = 0, lty = 2)

boxplot(E ~ CO$fHabitat, data = CO2)
abline(h = 0, lty = 2)

boxplot(E ~ CO$Country, data = CO2)
abline(h = 0, lty = 2)

#Not too bad. Though I still prefer to avoid a transformation
##############################################################






#Step 2: Is everything significant?
rownames(OUT1)[1:K] <- colnames(X)
OUT1

MyBUGSHist(out, 
           c(uNames("beta", K), "sigma", "sigma_Site"),
           PanelNames = MyNames)

#No...not everything is significant!
#But the terms that are not significant are levels from a factor!
#The credible intervals tell us whether they are different from each other.


#Options to figure out whether fHabitat is needed:
#1. Get DICs of model with and without fHabitat
    out$DIC  #56.01 with fHabitat versus 84.52 without habitat
#2. Get AICs of model with and without fHabitat
    #See the Bailey advanced code part
#3. Read the O’Hara and Sillanpaa (2009) paper and implement
#   one of their options
#4. Calculate the distribution of the differences between parameters

#Options 1/2 are relatively easy. Option 3 is the real way of doing it.
#Option 4 is very simple:

betas <- out$sims.list$beta
# These are the 6000 iterations for each parameter. We can 
# now calculate differences between parameters and 
# calculate a credible interval for the differences:
levels(CO$fHabitat)
#For example: FRS versus SHG
Dif <- betas[,4] - betas[,6]
hist(Dif)
#And get the credible interval
quantile(Dif, probs = c(0.025, 0.975))
#Is 0 in the interval? Yes...not significant.


#For example: APZ versus FRMZ
Dif <- betas[,1] - betas[,3]
quantile(Dif, probs = c(0.025, 0.975))
mean(Dif)

par(mfrow = c(3,1))
hist(betas[,1], xlim = c(0, 2))
hist(betas[,3], xlim = c(0, 2))
hist(Dif, xlim = c(0, 2))




# Repeat this for each combination of the fHabitat variables
# and make a table.
# Note...there is no need to correct for multiple testing
# There is no 5% chance that we make the wrong statement!


DoBayesianPosthoc <- TRUE
#Only run this code of the model is of the form
#intercept + LCC + fHabitat
if (DoBayesianPosthoc){

	betasHabitat <- out$sims.list$beta[,c(1,3:6)]
    A <- matrix(nrow = 5, ncol = 5)
    for (i in 1:5){
    ai <- betasHabitat[,i]   
	for (j in i:5){
        aj <- betasHabitat[,j]
        Difaiaj <- ai - aj
        qij <- quantile(Difaiaj, probs = c(0.025, 0.975))
        A[i,j] <- as.numeric(sign(qij[1] * qij[2]))
        A[j,i] <- A[i,j]
	}
}
A

#install.packages("corrplot")
library(corrplot)
col <- colorRampPalette(c("#BB4444", 
                         "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(A, 
         method = "shade", 
         diag = FALSE,
         shade.col = NA,
         tl.col = "black",
         tl.srt = 45,
         col = col(200))
#Blue means different from each other	
}
DoBayesianPosthoc <- FALSE


############################################################






########################################################
#Step 3: Explain the model
OUT1
#Write out the 5 equations....and sketch the fitted values

#E(LogG_ij) = 1.17 + 0.459 * LCC.std



##########################################################
#Final task....sketch the fitted values with CI
# We have a large number of betas
# We can calculate the fitted values for each MCMC
# iteration and use these to calculate credible intervals:

betas <- out$sims.list$beta

#Create covariate values for prediction
#But don't extrapolate!
tapply(CO$LCC.std, FUN = range, INDEX = CO$fHabitat)                     

LCC1 <- seq(-1.37, 3.533, length = 25)         
LCC2 <- seq(-1.19, 3.26, length = 25)         
LCC3 <- seq(-1.20, 1.68, length = 25)         
LCC4 <- seq(-1.18, -0.29, length = 25)         
LCC5 <- seq(-1.55, -0.54, length = 25)         

#Here is the covariate matrix
MyData <- data.frame(LCC.std = c(LCC1, LCC2, LCC3, LCC4, LCC5),
                     fHabitat = rep(levels(CO$fHabitat), each = 25))
              
              
              
#Here is some fancy code to do the same thing:
library(plyr)
MyData <- ddply(CO, .(fHabitat), summarize,
                LCC.std = seq(min(LCC.std), max(LCC.std),length = 10))
MyData

              
                     
                     
Xp     <- model.matrix(~ LCC.std + fHabitat, data = MyData)

#Now calculate 6000 fitted values
mu.mcmc <- Xp %*% t(betas)

#Calculate the posterior mean and 95% credible intervals
#based on the large number of lines
MyLinesStuff <- function(x){
   OUT <- matrix(nrow = nrow(x), ncol=4) 
	for(i in 1:nrow(x)){
	  xi <- x[i,]	
      OUT[i,3:4] <- quantile(xi, probs = c(0.025, 0.975))
      OUT[i,1] <- mean(xi)
      OUT[i,2] <- sd(xi)
	}
	colnames(OUT) <- c("mean", "se", "lb", "ub")
	OUT
}

L       <- MyLinesStuff(mu.mcmc)
MyData2 <- cbind(MyData,L)
MyData2

#Each row in L is an artificial covariate 
#Now plot the data and draw lines for the mean,
#and 2.5% and 97.5% values:


#Back transform the PD.std value

#           x - mean(x)
# x.std =  ------------- 
#            std(x)

MyData2$LCC <- MyData$LCC.std * sd(CO$LCC) + mean(CO$LCC)


#And use ggplot to make some nice graphs
library(ggplot2)
p <- ggplot()
p <- p + geom_point(data = CO, 
                    aes(y = LogG, x = LCC),
                    shape = 16, 
                    size = 3)
p <- p + xlab("LCC") + ylab("LogG")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData2, 
                   aes(x = LCC, y = mean), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData2, 
                     aes(x = LCC, 
                         ymax = ub, 
                         ymin = lb ),
                     alpha = 0.2)

p <- p + geom_hline(yintercept=log(3), lty = 2)                    

p <- p + facet_grid(. ~ fHabitat, scales = "fixed")
p

