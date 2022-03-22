#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#########################################################################
#Lilies exercise
#Set the working directory and import the data
setwd("/Users/Highstat/MaggieWindows/applicat/HighlandStatistics/webHighstatNew2_2010/CourseBaFrGLMM5/Data/")

LS <- read.table("Lilies.txt", 
                 header = TRUE)
##########################################################################
#Data from
#Beavers and lilies: selective herbivory and adaptive foraging
#behaviour
#Freshwater Biology (2014) 59, 224-232 doi:10.1111

#Based on pads collected from lakes where beavers are
#absent or where there is no evidence of grazing on
#N. alba, a strong positive relationship between
#the petiole diameter and midline distance of N. alba
#pads is to be expected!

names(LS)
# "Site"              "Mid_line_Diameter" "MD_um"
# "Petiole_Diameter"  "Beaver"


str(LS)
#'data.frame':   356 obs. of  5 variables:
# $ Site             : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Mid_line_Diameter: num  118.5 87.6 95.9 98.9 99.1 ...
# $ MD_um            : int  11851 8759 9594 9886 9908 7887 9684 8718 ...
# $ Petiole_Diameter : num  6.47 5.74 5.76 5.48 6.53 4.3 5.93 5.72 5.74 5.11 ...
# $ Beaver           : Factor w/ 2 levels "N","Y": 1 1 1 1 1 1 1 1 1 1 ...

######################################################################
# Aim:
# Model Mid_line_Diameter  as a function of:
#   -Petiole_Diameter
#   -Site effect

###################################################################
#Load packages and library files
#Change the path of this:
source(file = "/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV9.R")  
source(file = "/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstatV4.R")

library(lattice)
library(nlme)  #For lme
library(R2jags)
######################################################################



######################################################################
#Housekeeping
#Create variables with shorter names
LS$MD <- LS$Mid_line_Diameter
LS$PD <- LS$Petiole_Diameter

#Convert into factors
LS$fSite <- factor(LS$Site)
str(LS)

table(LS$fSite)
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
# 15 20 19 34 21 20 20 24 31 12 20 20 20 20 20 20 20
######################################################################




######################################################################
#Data exploration
## OUTLIERS?
par(mfrow = c(2, 1), mar = c(3, 4, 1, 1))
dotchart(LS$MD, main= "Mid_line_Diameter")
dotchart(LS$PD, main= "Petiole_Diameter")
#No!


######################################################################
# RELATIONSHIPS?
plot(x = LS$PD, 
     y = LS$MD,
     pch = 16, 
     cex = 0.5, 
     ylab = "Mid line Diameter (mm)",
     xlab = "Petiole Diameter (mm)")
#Perfect linear relationship!


#Site effect
boxplot(LS$MD ~ fSite,
        varwidth = TRUE, 
        data = LS,
        xlab = "Site",
        ylab = "Mid_line_Diameter")

#################################################################
#END OF DATA EXPLORATION
#################################################################




#Continue with the MCMC analysis
######################################################################
#Protocol for mixed modelling (3 steps).
#Based on Zuur et al. (2013).
#Beginner's Guide to GLM and GLMM with R

#1. Based on prior knowledge of the dependency
#   structure in the data, select the random
#   structure (i.e. the random effects) a priori.
#2. Fit the model and investigate which covariates
#   in the fixed part are significant or important.
#   Alternatively, omit model selection on the
#   covariates and adopt an information theoretic
#   approach (specify and compare 10 - 15 models).
#3. When the optimal model has been found, present
#   the numerical output and provide a graphic
#   representation of the model fit.

#In step 2 you may want to consider applying a limited amount of model
#selection...e.g. drop non-significant interactions.




##################################################################
#Step 1 of the protocol
#We recommend that you always try to compare frequentist and
#Bayesian models...just to exclude the chance that you make programming
#mistakes. But for more complicated models this is not possible as
#MCMC is quite often the only option


#Standardize each continuous covariate (A MUST!)
LS$PD.std <- (LS$PD - mean(LS$PD)) / sd(LS$PD)

#Model for comparing
M1 <- lme(MD ~ PD.std,
          random =~ 1 | fSite,
          data = LS, 
          method = "REML")

summary(M1)
##################################################



#Use MCMC to estimate the parameters
#First we will use lillies specific coding:
#1. Bundle data

#Random effects
#We need a vector with values 1 1 1 2 2 2 3 3 4 ...
Site <- as.numeric(as.factor(LS$fSite))
Site


# This is the model:
# MD_ij ~ N(mu_ij, sigma^2)
# E(MD_ij) = Fixed stuff + a_i 

# Fixed stuff = beta_1 +  PD.std * beta_2 
# i = 1, .., 17
# a_i ~ N(0, sigma_Site^2)
# We have 356 observations in total


# List with the response variable, covariates
# and site index
JAGS.data <- list(MD    = LS$MD,      #Response
                  PD.std = LS$PD.std, #Covariates
                  Site   = Site)      #Site index
JAGS.data



###################################################
#2. JAGS modelling code
#COPY AND PASTE in R FROM HERE.......
sink("liliesmixedmodel.txt")
cat("
model{
    #1A. Diffuse normal priors beta 
    # betas ~ N(0, 100^2)
    for (i in 1:2) { beta[i] ~ dnorm(0, 0.0001)}
    
    #1B. Diffuse uniform prior sigma 
    sigma ~ dunif(0, 20)
    tau  <- 1 / (sigma * sigma) 

    #1C. Priors random effects
    for (i in 1:17) { a[i] ~ dnorm(0, tau_Site) }
 
    #1D. Diffuse uniform prior for sigma Site
    tau_Site  <- 1 / (sigma_Site * sigma_Site)
    sigma_Site ~ dunif(0, 20)


    ######################## 
    #2. Likelihood (with lilies specific coding)
    for (i in 1:356) {
      MD[i]   ~ dnorm(mu[i], tau)
      eta[i] <- beta[1] + beta[2] * PD.std[i] #Covariates
      mu[i]  <- eta[i]  + a[Site[i]]          #Covariates + random effects
    }

    #Residuals    
    for (i in 1:356) {
      Res[i] <- MD[i] - mu[i] 
    }
}
",fill = TRUE)
sink()
#TO HERE...!!

#####################################



#3. Initial values & parameters to save
inits  <- function () {
  list(
    beta       = rnorm(2,  0, 0.1),
    sigma      = runif(1,  0, 20),
    a          = rnorm(17, 0, 0.1),
    sigma_Site = runif(1,  0, 20))  }

params <- c("beta", "sigma", "sigma_Site",
            "Res", "a",  
            "mu")


##############
#4. Start JAGS
J1 <- jags(data       = JAGS.data,
           inits      = inits,
           parameters = params,
           model      = "liliesmixedmodel.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

J2  <- update(J1, n.iter = 10000, n.thin = 10)
out <- J2$BUGSoutput



#5. Assess mixing and auto-correlation
MyNames <- c("Intercept", "PD", "sigma eps", "sigma site")
#Adjust this variable if extra parameters (e.g an ICC) are added!


MyBUGSChains(out, 
             c(uNames("beta", 2), "sigma", "sigma_Site"),
             PanelNames = MyNames)


MyBUGSACF(out, 
          c(uNames("beta", 2), "sigma", "sigma_Site"),
          PanelNames = MyNames)



#6. Present numerical output and posterior distribution
OUT1 <- MyBUGSOutput(out, 
                     c(uNames("beta", 2), "sigma", "sigma_Site"),
                     VarNames = MyNames)
print(OUT1, digits = 5)

MyBUGSHist(out, 
           c(uNames("beta", 2), "sigma", "sigma_Site"),
           PanelNames = MyNames)


#Compare Bayesian and frequentist results
OUT1
summary(M1)
######################################################





######################################################
#Model validation

#Extract residuals and fitted values from JAGS
E1.mcmc <- out$mean$Res  #Posterior mean of the residuals
F1.mcmc <- out$mean$mu   #This is with random effects

# So...each value represents the mean of a large number
# of MCMC iterations



#1.1 - Plot residuals versus fitted values
#1.2 - Plot residuals versus each covariate in the model
#1.3 - Plot residuals versus each covariate NOT in the model
#1.4 - Plot residuals versus time (if relevant)
#1.5 - Plot residuals versus spatial coordinates


#1.1
par(mfrow = c(1,1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1.mcmc, 
     y = E1.mcmc,
     xlab = "Posterior mean values",
     ylab = "Posterior mean residuals")
abline(h = 0, lty = 2)


#1.2
par(mfrow = c(1,1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
boxplot(E1.mcmc ~ fSite, 
        data = LS,
        ylab = "Posterior mean residuals",
        xlab = "Site")
abline(h = 0, lty = 2)


par(mfrow = c(1,1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = LS$PD, 
     y = E1.mcmc,
     xlab = "PD",
     ylab = "Posterior mean residuals")
abline(h = 0, lty = 2)



#1.3: Not relevant here

#1.4: Not relevant here .... is it?

#1.5 : Not relevant here

#Results model validation:
#Looks ok (?)...no structure in residuals...go to step 2
#########################################################




########################################
#Step 2: Do we need all covariates? Is everything significant?
print(OUT1, digits = 5)

MyNames <- c("Intercept", "PD", "sigma eps", "sigma site")

MyBUGSHist(out, 
           c(uNames("beta", 2), "sigma", "sigma_Site"),
           PanelNames = MyNames)


#Clearly PD is significant.
#################################################################






#################################################################
# Model interpretation

# Task:  Explain what the model tells you.
OUT1


#Sketch the fitted values with CI
range(LS$PD.std)   #-2.337537  2.479650
MyData <- data.frame(PD.std = seq(-2.34, 2.48, length = 25))
Xp     <- model.matrix(~ PD.std, data = MyData)
dim(Xp)


# We have a large number of MCMC values of the betas.
# We can calculate the fitted values for each MCMC iteration  
# and use these to calculate credible intervals:

#Get the MCMC betas
Beta.mcmc <- out$sims.list$beta 
dim(Beta.mcmc)


# Calculate the fitted values for each MCMC iteration
# These are fitted values for a typical site 
# (typical = average = 0)

mu.mcmc <- Xp %*% t(Beta.mcmc)  #25 by 2 times 2 by 3000   
dim(mu.mcmc)

# Note that this mu.mcmc contains predicted values
# for the 25 artifical covariate values.
# Now what?

#We could plot these....
par(mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = LS$PD.std, 
     y = LS$MD, 
     xlab = "Standarized Petiole_Diameter",
     ylab = "MD", 
     pch = 16)
for (i in 1:100){      
  lines(x = MyData$PD.std, 
        y = mu.mcmc[,i], 
        col = i, 
        lty = 1)
}
# We only plotted fitted values for the first 100 
# MCMC iterations..

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


par(mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = LS$PD.std, 
     y = LS$MD, 
     xlab = "Standarized Petiole_Diameter",
     ylab = "MD", 
     pch = 16)

lines(MyData$PD.std, L[,"mean"], col = 2, lwd = 2)
lines(MyData$PD.std, L[,"lo"],   col = 2, lwd = 2, lty = 2)
lines(MyData$PD.std, L[,"up"],   col = 2, lwd = 2, lty = 2)


# If you would like to make this plot in ggplot2, then 
# it is easier to have L and MyData combined:

MyData2 <- cbind(MyData,L)
MyData2

# And perhaps you would like to have the original PD
# values and not the standardized values.
# This is how we standardized the data:

#           x - mean(x)
# x.std =  ------------- 
#              std(x)

# To back-standardize, use:
MyData2$PD <- MyData$PD.std * sd(LS$PD) + mean(LS$PD)

# The rest is some fancy ggplot2 coding:
library(ggplot2)
p <- ggplot()
p <- p + geom_point(data = LS, 
                    aes(y = MD, x = PD),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Petiole Diameter (mm)") + ylab("Mid line Diameter (mm)")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_line(data = MyData2, 
                   aes(x = PD, y = mean), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData2, 
                     aes(x = PD, 
                         ymax = up, 
                         ymin = lo ),
                     alpha = 0.4)
p
########################################
#End of visualisation
      
      
      
      
#Task: Question from a referee: I would like to see
#      the posterior distribution of the ICC
#      To answer this question, adjust the JAGS code 
#      above.




      
      
      
#####################################################
#Bayesian post-hoc test on the random effects
#Are the random effects different from each other?
out$sims.matrix[,"a[1]"]  

# Get all MCMC values for the random
# intercepts at site 1 and site 2:
a1 <- out$sims.matrix[,"a[1]"]   #iterations for a1
a2 <- out$sims.matrix[,"a[2]"]   #iterations for a2

#Calcuate the difference and show the posterior
#distribution of the difference
Difa1a2 <- a1 - a2               #Difference
hist(Difa1a2, xlim = c(-40, 5))  #Histogram
abline(v = 0, col = 2, lty = 2, lwd = 2)
q12 <- quantile(Difa1a2, probs = c(0.025, 0.975))
q12


#Do the same for other sites
A <- matrix(nrow = 17, ncol = 17)
for (i in 1:17){
    ai <- out$sims.matrix[,paste("a[", i, "]", sep = "")]   
	for (j in i:17){
        aj <- out$sims.matrix[,paste("a[", j, "]", sep = "")]   
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


#Limitation of this solution file:
#Everything is lillie-specific.

         