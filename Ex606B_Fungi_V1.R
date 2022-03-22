
#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
# Fungi exercise
# These data are taken from:
# Elucidating the nutritional dynamics of fungi using 
# stable isotopes
# Jordan R. Mayor, Edward A. G.Schuur and Terry W. Henkel
# Ecology Letters, (2009) 12: 171-183


#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/NewDataSets/lme_Fungi/")
Fungi <- read.csv("FungiV5.csv", 
                   header = TRUE,
                   dec = ".")
names(Fungi)
 # [1] "Study"    "Species"  "Site"     "Temp"     "Rain"     "Type"    
 # [7] "d13C"     "d15N"     "Category" "STND_13c" "STND_15n" "Lat"     
# [13] "Long"    

str(Fungi)

# Fungi function at two fundamental biogeochemical 
# interfaces between soil and plants. 
#   -Decomposer fungi mineralize organic carbon compounds 
#    in detritus and liberate mineral nutrients in the process.
#   -Mycorrhiza-forming fungi, as mutualistic root extensions, 
#    enhance mineral, and perhaps organic, nutrient uptake in 
#    exchange for plant photosynthate

# Fungi are divided into saprotrophic (SAP) and ectomycorrhizal (ECM) 
# functional groups

# Mycorrhizal and saprotrophic (SAP) fungi are essential 
# to terrestrial element cycling due to their uptake of 
# mineral nutrients and decomposition of detritus. 

# d15N and d13C of fungi provide time integrated biogeochemical# information regarding the acquisition, transformation,# and export of C and N by fungi.

# In order to determine if fungal d15N or d13C# patterns across ecosystems are similar to plants and 
# soils, we assessed the explanatory capacity of:
#   -mean annual temperature (Temp...called MAT in paper),
#   -mean annual precipitation (Rain..called MAP in paper), and 
#   -latitude (Lat). 

#Methods:
# 1. The athors compiled d15N and d13C values from 
#    one novel and ten published data sets. Compiled 
#    data included 913 d15N and 813 d13C values from 
#    collector categorized ECM or SAP fungi, and 27 
#    fungi of unknown ecological role, together comprising 
#    148 genera.
#2.  32 study sites
###################################



###########################################################
#Load packages and library files
library(lattice)  #Needed for multi-panel graphs
library(nlme)
library(mgcv)
library(ggplot2)
library(R2jags)

#source the file: HighstatLibV7.R
source("/Users/Highstat/MaggieWindows/applicat/HighlandStatistics/webHighstatNew2_2010/CourseMCMCGLMM/RCode/HighstatLibV7.R")

source(file = "/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstatV3.R")

##################################################################


# Task:
# Model d13C and d15N as a function of:
# Temp, Rain, Lat, Type (ECM, SAP, Unknown) + relevant interactions
# Random effect site 
###################################################




###################################################
#Housekeeping
Fungi$fSite    <- factor(Fungi$Site)

###################################################



###################################################
#Data exploration
#Outliers
MyVar <- c("Temp", "Rain", "Lat", "Long", "d15N", "d13C")
Mydotplot(Fungi[,MyVar])
# One small d13C value
# Lots of observations with the same covariate values!

#Are the data balanced?
table(Fungi$fSite)
sort(table(Fungi$Species))
table(Fungi$Study)
table(Fungi$Type)


#Collinearity
boxplot(Temp ~ Type, 
        data = Fungi, 
        xlab = "Type", 
        ylab = "Temp")
#Small collinearity?

boxplot(Rain ~ Type, 
        data = Fungi, 
        xlab = "Type", 
        ylab = "Rain")
#Minor collinearity?

#Site will be collinear with temp and rain.
#And latitude


MyVar <- c("d13C", "d15N", "Temp", "Rain","Lat")
Mypairs(Fungi[, MyVar])
#Serious collinearity! Pick either Temp, rain or lat!
#Pick Rain?

p <- ggplot()
p <- p + geom_point(data = Fungi, 
                    aes(y = d15N, x = Rain),
                    shape = 16, 
                    size = 3)
p <- p + xlab("Rain") + ylab("d15N")
p <- p + theme(text = element_text(size=15))
p <- p + geom_smooth(data = Fungi, method = "lm", 
                   aes(x = Rain, y = d15N), 
                   colour = "black")                     
p <- p + facet_grid(. ~ Type, scales = "fixed")
p

#There is probably an interaction. But the quality of the
#rain data in the unknown group is not that nice!
#Should we analyse the data without the unknowns?


####################
#Relationships
boxplot(d15N ~ fSite, 
        data = Fungi, 
        xlab = "Site")
#Site effect!
        
boxplot(d15N ~ Study, 
        data = Fungi, 
        xlab = "Study") 
#Location effects?

boxplot(d15N ~ Species, 
        data = Fungi, 
        xlab = "Species")
#Not many observations per species


par(mar = c(5,6,2,2), cex.lab = 1.5)
plot(x = Fungi$Rain, 
     y = Fungi$d15N,
     pch = 16,
     cex = 1,
     xlab = "Rain",
     ylab = "d15N")

par(mar = c(5,6,2,2), cex.lab = 1.5)
plot(x = Fungi$Temp, 
     y = Fungi$d15N,
     pch = 16,
     cex = 1,
     xlab = "Temperature",
     ylab = "d15N")


par(mar = c(5,6,2,2), cex.lab = 1.5)
plot(x = Fungi$Lat, 
     y = Fungi$d15N,
     pch = 16,
     cex = 1,
     xlab = "Latitude",
     ylab = "d15N")
     
     
     
#Dump NAs
MyVar <- c("d13C", "d15N", "Temp", "Rain", 
           "Lat", "Long", "fSite", "Type",
           "Study", "Species")
Fungi2 <- na.exclude(Fungi[, MyVar])
dim(Fungi)
dim(Fungi2)


table(Fungi2$fSite, Fungi2$Type)
#it would make sense to do only the ECM and SAP data!
#############


#Data exploration indicates that:
# There are is potentially 1 outlier
# There is serious collinearity
# There could be an interaction
# Should we drop the Unknowns?


##########################################################



#Start of analysis. Bayesian analysis
#Drop the unknowns
Fungi3 <- Fungi2[Fungi2$Type != "Unknown", ]
table(Fungi3$Type)
Fungi3$Type <- factor(Fungi3$Type)

Fungi3$Rain.std <- as.numeric(scale(Fungi3$Rain)) 
#Reference model
M0 <- lme(d13C ~ Rain.std *  Type,
          random =~ 1 | fSite,
          data = Fungi3, 
          method = "REML")
summary(M0)
#By design of the study we use site as random effect.
      
      
      
#Use MCMC to estimate the parameters
#First we will use lillies specific coding:
#1. Bundle data

X <- model.matrix(~ Rain.std *  Type, data = Fungi3)
K <- ncol(X)  #Number of betas
head(X)

#Random effects
#We need a vector with values 1 1 1 2 2 2 3 3 4 ...
Site <- as.numeric(factor(Fungi3$fSite))
Site
Nre <- length(unique(Site))


#TASK:  Write down the model that we are going to fit
#


# List with the response variable, covariates
# and site index
win.data <- list(Y    = Fungi3$d13C, #Response
                 X    = X,     #Covariate
                 K    = K,     #Number of betas
                 N    = nrow(Fungi3), #Sample size
                 Site = Site,  #Random effect ID
                 Nre  = Nre)   #Number of random effects
win.data



###################################################
#2. JAGS modelling code
sink("fungimm.txt")
cat("
model{
    #1A. Diffuse normal priors beta 
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}
    
    #1B. Diffuse uniform prior sigma
    sigma ~ dunif(0, 20)
    tau  <- 1 / (sigma * sigma) 

    #1C. Priors random intercepts and random slopes
    for (i in 1:Nre) { a[i] ~ dnorm(0, tau1_Site) }

    #1D. Diffuse uniform prior for sigma1_Site
    #    This one is for the random intercept
    tau1_Site  <- 1 / (sigma1_Site * sigma1_Site)
    sigma1_Site ~ dunif(0, 100)


    ######################## 
    #2. Likelihood 
    for (i in 1:N) {
      Y[i]    ~ dnorm(mu[i], tau)
      eta[i] <- inprod(beta[], X[i,]) #Covariates
      mu[i]  <- eta[i]  + a[Site[i]]  #Covariates + random intercepts 
    }

    #Residuals    
    for (i in 1:N) {
      Res[i]    <- Y[i] - mu[i] 
    }
}
",fill = TRUE)
sink()
#####################################



#3. Initial values & parameters to save
inits  <- function () {
  list(
    beta       = rnorm(K,  0, 0.1),  #betas
    sigma      = runif(1,  0, 20),   #sigma eps
    a          = rnorm(Nre, 0, 0.1),  #random intercepts
    sigma1_Site = runif(1,  0, 100))  } #sigma random intercepts

params <- c("beta", "sigma", "sigma1_Site",
            "Res", "a",  
            "mu")


##############
#4. Start JAGS
G1 <- jags(data       = win.data,
           inits      = inits,
           parameters = params,
           model      = "fungimm.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

G2   <- update(G1, n.iter = 10000, n.thin = 10)
outG <- G2$BUGSoutput



#5. Assess mixing and auto-correlation
MyNames <- c("Intercept", "Rain", "Type", "Interaction",
             "sigma eps", "sigma1 site")
#Adjust this variable if extra parameters (e.g an ICC) are added!


MyBUGSChains(outG, 
             c(uNames("beta", K), "sigma", "sigma1_Site"),
             PanelNames = MyNames)


MyBUGSACF(outG, 
          c(uNames("beta", K), "sigma", "sigma1_Site"),
          PanelNames = MyNames)



#6. Present numerical output and posterior distribution
OUT3 <- MyBUGSOutput(outG, 
                     c(uNames("beta", K), "sigma", 
                       "sigma1_Site"),
                     VarNames = MyNames)
print(OUT3, digits = 5)

MyBUGSHist(outG, 
           c(uNames("beta", K), "sigma", "sigma1_Site"),
           PanelNames = MyNames)


#Compare Bayesian and frequentist results
OUT3
summary(M1)
######################################################



#Task:  Write down the fitted model




######################################################
#Model validation

#Extract residuals and fitted values from JAGS
E1.mcmc <- outG$mean$Res  #Posterior mean of the residuals
F1.mcmc <- outG$mean$mu   #This is with random effects

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
        data = Fungi3,
        ylab = "Posterior mean residuals",
        xlab = "Site")
abline(h = 0, lty = 2)


par(mfrow = c(1,1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = Fungi3$Temp, 
     y = E1.mcmc,
     xlab = "Temp",
     ylab = "Posterior mean residuals")
abline(h = 0, lty = 2)


#TASK:
#Make the other model validation graphs



#1.3: Not relevant here

#1.4: Not relevant here

#1.5 : Not relevant here

#Results model validation:
#Looks ok (?)...no structure in residuals...go to step 2
#########################################################




########################################
#Step 2: Do we need all covariates? Is everything significant?
print(OUT3, digits = 5)

MyBUGSHist(outG, 
           c(uNames("beta", K), "sigma", "sigma1_Site"),
           PanelNames = MyNames)

#################################################################






#################################################################
# Model interpretation

# Task:  Explain what the model tells you.
OUT3






#Sketch the fitted values with CI
library(plyr)
MyData <- ddply(Fungi3, .(Type), summarize,
                Rain.std = seq(min(Rain.std), max(Rain.std),length = 25))
MyData
X  <- model.matrix(~ Rain.std *  Type, data = MyData)


# We have a large number of MCMC values of the betas.
# We can calculate the fitted values for each MCMC iteration  
# and use these to calculate credible intervals:

#Get the MCMC betas
Beta.mcmc <- outG$sims.list$beta 
dim(Beta.mcmc)


# Calculate the fitted values for each MCMC iteration
# These are fitted values for a typical site 
# (typical = average = 0)

mu.mcmc <- X %*% t(Beta.mcmc)  #25 by 4 times 4 by 3000   
dim(mu.mcmc)

# Note that this mu.mcmc contains predicted values
# for the 25 artifical covariate values.

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

# And perhaps you would like to have the original PD
# values and not the standardized values.
# This is how we standardized the data:

#           x - mean(x)
# x.std =  ------------- 
#              std(x)

# To back-standardize, use:
MyData2$Rain <- MyData$Rain.std * sd(Fungi3$Rain) + mean(Fungi3$Rain)

# The rest is some fancy ggplot2 coding:
p <- ggplot()
p <- p + geom_point(data = Fungi3, 
                    aes(y = d13C, x = Rain),
                    shape = 1, 
                    size = 1)
p <- p + xlab("rain") + ylab("d13C")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_line(data = MyData2, 
                   aes(x = Rain, y = mean), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData2, 
                     aes(x = Rain, 
                         ymax = up, 
                         ymin = lo ),
                     alpha = 0.4)
p <- p + facet_grid(. ~ Type, scales = "fixed")
p


#Or:
p <- ggplot()
p <- p + geom_point(data = Fungi3, 
                    aes(y = d13C, 
                        x = Rain,
                        colour = Type),
                    shape = 1, 
                    size = 1)
p <- p + xlab("rain") + ylab("d13C")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_line(data = MyData2, 
                   aes(x = Rain, 
                       y = mean, 
                       group = Type,
                       colour = Type))

p <- p + geom_ribbon(data = MyData2, 
                     aes(x = Rain, 
                         ymax = up, 
                         ymin = lo,
                         group = Type,
                         colour = Type,
                         fill = Type ),
                     alpha = 0.4)
p

#Task: Explain what it all means