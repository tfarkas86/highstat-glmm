#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
#Cowbird exercise
#################################################################
#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/FilesMCMC_GLMM/CDMCMCGLMGLMMCourse/Data")
CB <- read.table(file = "CowbirdV3Book.txt",
                 header = TRUE)
dim(CB)
names(CB)
str(CB)
#################################################################

#Task: 
#Some bird species place their eggs solely in the nests 
#of other bird species; this is called brood parasitism. 
#Brood parasitism has led to an evolutionary battle between 
#host and parasite; the parasite tries to fool the host 
#(e.g., mimic the host eggs), and the host tries to outsmart 
#the parasite with advanced defensive mechanisms (e.g., nest location). 
#Patten et al. (2011) used an extensive dataset and looked for 
#cues in brood parasitism nest selection by the brood parasitic 
#brown-headed Cowbird (Molothrus ater).

#The dataset consists of a large number of observations (2,825) 
#on absence and presence of brood parasitism by cowbirds. 
#Patten et al. (2011) used logistic regression to investigate 
#whether the probability of brood parasitism is affected by 
#covariates such as perch proximity, nest height, livestock 
#proximity, habitat density, nest exposure, laying date of 
#the first egg, host clutch size, and host species, among other 
#variables. They found egg laying date to have a significant effect.
#The question we focus on in this chapter builds on the work 
#by Patten and colleagues: 

#Is the effect of egg laying date on the probability of brood 
#parasitism by cowbirds the same for all host species, or does this 
#effect differ per host species? 
#And if the effect differs per species, how does it differ from 
#species to species?

###################################################################
#Load packages and library files
library(lattice)  #Needed for multi-panel graphs
library(lme4)
source(file="HighstatLibV9.R")  

source("MCMCSupportHighstatV4.R")

library(R2jags)

########################################################




########################################################
#Data exploration

#First we dump a few NAs
CB1 <- na.exclude(CB)
dim(CB1)

#Outliers
dotchart(CB1$Firstegg)
#That is ok

table(CB1$Cowbird)
table(CB1$Plot)
table(CB1$Year)
#That is all ok

table(CB1$Species)
#This is not ok!! If we want to include interactions, then
#we should dumpt all species < 25 or 30-ish obervations

#These data were used for a GAMM chapter, and we included
#a Firstegg x Species smoother interaction. For this 
#you want to have at least 50-ish or more number
#of observations per species.

#This is whaty we used in the 2014 book
CB2 <- subset(CB1, Species == "DICK" |
                   Species == "EAME" |
                   Species == "GRSP" |
                   Species == "RWBL")
CB2$Species <- factor(CB2$Species)
#But for GLMM techniques you could as well
#add COGR, EAPH and OROR
#Using only 4 species keeps computing time
#and coding of GAMM relatively low and simple.


#Check whether this is still ok:
table(CB2$Cowbird)
table(CB2$Plot)
table(CB2$Year)
#Yes...that is all ok

     
#Collinearity
boxplot(Firstegg ~ Year, data = CB2) 
boxplot(Firstegg ~ Species, data = CB2) 
boxplot(Firstegg ~ Year, data = CB2) 
#Some minor collinearity going on!

    
#Relationships    
plot(x = CB2$Firstegg,
     y = CB2$Cowbird)
#That is not very useful.
#Try this:
boxplot(Firstegg ~ Cowbird, data = CB2)     

with(CB2, table(Species, Cowbird))
with(CB2, table(Year, Cowbird))


           
#Balanced data?     
tapply(CB2$Firstegg, INDEX = CB2$Species, FUN = range)
# $DICK
# [1]  45 133
# $EAME
# [1]  19 114
# $GRSP
# [1]  34 131
# $RWBL
# [1]  26 116

#That could be trouble!
#Since we are interested in the interaction
#between species and Firstegg, it may be an option
#to focus the analysis on on the Firstegg data between
#days 46 and 116
#CB3 <- subset(CB2, Firstegg >= 46 & Firstegg <= 116)

#Alterntively..use the data as it is...and be careful 
#with the interpretation. That's what we will do here
CB3 <- CB2



#If you have a slow computer:
dim(CB3)
Computer <- "NotSlow"
if (Computer == "Slow") {
	set.seed(12345)
	CB3.backup <- CB3  #Make backup
	I <- sample(1:nrow(CB3), size = 300, replace = FALSE)
	CB.temp <- CB3[I,] 
    CB3     <- CB.temp  #Overwrite existing file
   }
Computer <- "Dont run again"
dim(CB3)





################################################
#Some housekeeping
CB3$fSpecies <- factor(CB3$Species)
CB3$fYear    <- factor(CB3$Year)
CB3$fPlot    <- factor(CB3$Plot)





#######################################################
#By design of the study use fPlot as random effect
#Step 1: Apply a binomial GLMM using MCMC in JAGS

#What will you learn in this excercise?
#1. How to do a Bernoulli GLMM in JAGS
#2. One of the major benefits of MCMC:
#   Everything you calculate comes with a whole
#   distribution...e.g. differences between fitted 
#   values.


#Standardize the covariate Firstegg
CB3$Firstegg.std <- (CB3$Firstegg - mean(CB3$Firstegg)) / sd(CB3$Firstegg)

#Reference model:
#detach(package:nlme)
M1 <- glmer(Cowbird ~ Firstegg.std * fSpecies + fYear + (1|fPlot),
            data = CB3, 
            family = binomial)
summary(M1)



####################################
##Step 1: Fit the model.
#Use MCMC to estimate the parameters                     
#1. Bundle data
X <- model.matrix(~ Firstegg.std * fSpecies + fYear, data = CB3)                      
K <- ncol(X)  #Number of columns
head(X)

#Random effects:
Plot <- as.numeric(as.factor(CB3$fPlot))
Plot
Nre <- length(unique(CB3$fPlot))


#The code below is copied an pasted from the owls
#Nest was changed into Plot
win.data <- list(Y     = CB3$Cowbird, 
                 X     = X,
                 N     = nrow(CB3),
                 K     = ncol(X),
                 Plot  = Plot,
                 Nre   = Nre)
win.data

###################################################
# 2. JAGS modelling code
sink("GLMMCBJAGS.txt")
cat("
model{
    #1A. Priors beta and sigma
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}

    #1B. Priors random effects and sigma_Plot
    tau_Plot <- 1 / (sigma_Plot * sigma_Plot)
    
    for (i in 1:Nre) { a[i] ~ dnorm(0, tau_Plot)}
    
    num           ~ dnorm(0, 0.0016) # 1 / 25^2      
    denom         ~ dnorm(0, 1)               
    sigma_Plot   <- abs(num / denom)          

    #2. Likelihood
    for (i in 1:N) {
      Y[i]          ~ dbern(Pi[i])  
      logit(Pi[i]) <- eta[i] 
      eta[i]       <- inprod(beta[], X[i,]) + a[Plot[i]] 
  
      #Saving computer time
      ##3. Discrepancy measures 
      #Exp[i] <- Pi[i] 
      #Var[i] <- Pi[i] * (1 - Pi[i])
      #E[i]   <- (Y[i]  - Exp[i]) / sqrt(Var[i])    
        }          
}
",fill = TRUE)
sink()
#####################################
#Tools for assessing model fit are discussed in 
#Roos M, Held L (2011) Sensitivity analysis in Bayesian 
#generalized linear mixed models for binary data. 
#Bayesian Analysis 6: 259-278



#3. Initial values & parameters to save
inits  <- function () {
  list(
    beta  = rnorm(ncol(X), 0, 0.01),
    a     = rnorm(Nre, 0, 0.01),
    num   = rnorm(1, 0, 25),   
    denom = rnorm(1, 0, 1)  )  }

#This now becomes an issue for this data size.
#Let's take the minimal output
params <- c("beta", "a", "sigma_Plot")



#Don't forget to change the file name!
#4. Start JAGS
J1 <- jags(data       = win.data,
           inits      = inits,
           parameters = params,
           model      = "GLMMCBJAGS.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 1000,
           n.iter     = 3000)

#J2  <- update(J1, n.iter = 1000, n.thin = 1)  
J2  <- update(J1, n.iter = 10000, n.thin = 10)  

# J3 <- update(J2, n.iter = 10000, n.thin = 10)
# J4 <- update(J3, n.iter = 10000, n.thin = 10)


save(J2, file = "CowBirds_J2.RData")
#To reload results from an earlier session, use:
#load(file = "CowBirds_J2.RData")
#recompile(J2)
#load.module("dic")


out <- J1$BUGSoutput

print(J1, digits = 3)  #Takes a few minutes

#Write down DIC at this point


#5. Assess mixing
MyNames <- c(colnames(X), "Sigma plot")

MyBUGSChains(out, 
             c(uNames("beta", K), "sigma_Plot"), 
             PanelNames  = MyNames)
MyBUGSACF(out, 
          c(uNames("beta", K), "sigma_Plot"),
          PanelNames  = MyNames)


#6. Present output
K <- ncol(X)
OUT1 <- MyBUGSOutput(out, 
                    c(uNames("beta", K), "sigma_Plot"),
                    VarNames  = MyNames)
print(OUT1, digits =5)

MyBUGSHist(out, 
           c(uNames("beta", K), "sigma_Plot"),
           PanelNames  = MyNames)
#Some parameters are not significant!


OUT1
summary(M1)  #That is all good!
######################################################



#Now let's do something nice and show the benefits of MCMC.
#Let's calculate the fitted values for each species for 1 year,
#and subtract the fitted values from each other. 
#And also obtain a 95% credible interval for the difference. 
#Plot these differences, and the 95% credible interval versus 
#Firstegg, and see when it is different from 0
#So in MCMC language:

#1. Get the 6,000 betas
#2. Define a grid of covariate values
#3. Calculate 6,000 times the predicted values on this grid.
#4. Subtract the 6,000 predicted values for each species
#5. Caculate the 95% credible interval of the differences.
#6. Inspect whether 0 is in the 95% credible interval.

#In a frequentist setting it would be very difficult to calculate
#the distribution of the difference of two probabilities.
#In MCMC it is easy.

#1. Get the 6,000 betas
Beta.mcmc <- out$sims.list$beta  
dim(Beta.mcmc)

#2. Define a grid of covariate values
range(CB3$Firstegg.std)
MyData <- expand.grid(Firstegg.std = seq(-2.31, 2.74, length = 25),
                     fSpecies     = levels(CB3$fSpecies),
                     fYear        = factor(1992, 
                                     levels = c("1992", "1993", 
                                        "1994", "1995", "1996"))
                     )
Xp     <- model.matrix(~  Firstegg.std * fSpecies + fYear, 
                       data = MyData)

#3. Calculate 6,000 times the predicted values on this grid.
MyData$eta.mcmc <- Xp %*% t(Beta.mcmc)
dim(Xp)
dim(Beta.mcmc)
dim(MyData$eta.mcmc)

MyData$Pi <- exp(MyData$eta.mcmc) / (1 + exp(MyData$eta.mcmc))
dim(MyData$Pi)


#4. Subtract the 6,000 predicted values for each species
#"DICK" "EAME" "GRSP" "RWBL"

Pi.DICK <- MyData[MyData$fSpecies == "DICK", "Pi"]
Pi.EAME <- MyData[MyData$fSpecies == "EAME", "Pi"]
Pi.GRSP <- MyData[MyData$fSpecies == "GRSP", "Pi"]
Pi.RWBL <- MyData[MyData$fSpecies == "RWBL", "Pi"]

dim(Pi.DICK)
dim(Pi.EAME)
dim(Pi.GRSP)
dim(Pi.RWBL)

#So...these are the 6,000 fitted values at the 25 artificial Firstegg values
#Calculate differences:
Dif.GRSP.RWBL <- Pi.GRSP - Pi.RWBL

#5. Caculate the 95% credible interval of the differences.
dim(Dif.GRSP.RWBL)


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
L <- MyLinesStuff(Dif.GRSP.RWBL) 
L  #Posterior mean, se and 95% CI for each of these covariate values




#6. Inspect whether 0 is in the 95% credible interval.
FirstEgg.25 <- seq(-2.31, 2.74, length = 25)
plot(x = FirstEgg.25 * sd(CB3$Firstegg) + mean(CB3$Firstegg),
     y = L[,"mean"],
     lty = 1,
     type = "l",
     xlab = "Firstegg",
     ylab = "Difference between Prob GRSP and RWBL",
     ylim = c(-0.2, 0.2))
abline(h = 0, lty = 2)
lines(x = FirstEgg.25 * sd(CB3$Firstegg) + mean(CB3$Firstegg),
     y = L[,"lb"],
     lty = 2)
lines(x = FirstEgg.25 * sd(CB3$Firstegg) + mean(CB3$Firstegg),
     y = L[,"ub"],
     lty = 2)
  
#So..there are differences between these two species
#between days 45-ish and 80-ish. We can even get
#CIs for these days using MCMC!!
     
#Home work:
# 1. Do this for each species combination
# 2. Limit the Firstegg gradient to values
#    for which both species are sampled!

    













