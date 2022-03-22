
#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#######################################################################


#These data were taken from:
#Large Roads Reduce Bat Activity across Multiple Species
#Justin Kitzes, Adina Merenlender
#PLOS ONE | www.plosone.org 1 May 2014 | Volume 9 | Issue 5 | e96341

# The data:
# To examine the effects of large roads on bat populations,
# Kitzes and Merenlender (2014), used acoustic recorders to
# survey bat activity along ten 300 m transects bordering three
# large highways in northern California.

# Nightly counts of the hoary bat passes (Lasiurus cinereus),
# were analyzed to determine the relationship between bat activity
# and distance from a road. Sampling was conducted for a total
# of 34 nights in August and September of 2010 and 2011, with
# 5-6 points sampled on each night.

# Explanatory variables included three main variables of interest
# that could potentially influence total bat activity:
#  - distance from road (dist road, three levels: 0 m, 100 m, and 300 m),
#  - presence of a light within 100 m of a sampling point
#    (light 100, two levels: False, True), and 
#  - daily maximum temperature (temp max, continuous).

# In addition to these three main variables, four additional 
# variables were also recorded:
#  -site (DOED, HAYW, SAPA),
#  -year (2010, 2011), 
#  -night (lots of nights), and 
#  -transect (lots of transects) 

# Aim: Model LACI counts as a function of these covariates 
#      and figure out which are important
#######################################################################



#######################################################################
#Import the data
#For a Mac:
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data")
Bats <- read.table("BatsV2.txt", header = TRUE)

names(Bats)
str(Bats)

# [1] "Day"      "Month"    "Year"     "Site"     "Transect" "Road"    
# [7] "DistRoad" "RoadSide" "Light100" "TempMean" "TempMax"  "TempMin" 
#[13] "WindMean" "WindMax"  "ANPA"     "EPFU"     "MYEV"     "MYTH"    
#[19] "LANO"     "LACI"     "MYCA"     "LABL"     "PAHE"     "MYLU"    
#[25] "MYVO"     "TABR"     "MYYU"     "COTO"     "EUPE"     "TotAbund"
###########################################################################




###########################################################################
#Load support files and packages
source(file = "HighstatLibV9.R")
library(lme4)
library(lattice)
library(R2jags)
source("/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstatV3.R")


###########################################################################





###########################################################################
#Housekeeping
Bats$fDistRoad <- factor(Bats$DistRoad)
Bats$fLight100 <- factor(Bats$Light100)
Bats$fSite     <- factor(Bats$Site)
Bats$fYear     <- factor(Bats$Year)
Bats$fTransect <- factor(Bats$Transect)

###########################################################################



###########################################################################
#Data exploration

#Outliers
MyVar <- c("LACI", "Day", "Month", "Year",
           "DistRoad", "TempMean", "TempMax",  
           "TempMin", "WindMean", "WindMax")
            
Mydotplot(Bats[,MyVar])

#LACI: That will most likely be a Poisson
#Day: If we use this as random effect, then we need to use the 
#     month-day combination to define a unique day
#


###########################################################################
#Zero inflation
sum(Bats$LACI == 0) / nrow(Bats)
#30%


###########################################################################
# Are the categorical variables balanced?
# Define day as 'Day in Year'
Bats$Date <- paste(Bats$Day, Bats$Month, Bats$Year, sep = "/")
Bats$DayInYear <- strptime(Bats$Date, "%d/%m/%Y")$yday + 1
Bats$DayInYear
cbind(Bats$DayInYear, Bats$Day, Bats$Month, Bats$Year)

table(Bats$Date)
#Use this as a random effect to capture correlation within the same night
Bats$fDate <- factor(Bats$Date)


table(Bats$Year)
table(Bats$Site)
table(Bats$DistRoad)
table(Bats$DistRoad, Bats$Site) #For 2-way interaction
table(Bats$Transect) #Use this as random effect
# All ok-ish!
##########################################################################



##########################################################################
#Collinearity
MyVar <- c("DayInYear", "Month", "Year",
           "DistRoad", "TempMean", "TempMax",  
           "TempMin", "WindMean", "WindMax")
            
Mypairs(Bats[,MyVar])

#Use TempMax and WindMax
MyVar <- c("DayInYear", "Month", "Year",
           "DistRoad", "TempMax", 
           "WindMax")
            
Mypairs(Bats[,MyVar])

#Collinearity between continuous variable and factor
boxplot(TempMax ~ fYear, data = Bats)
boxplot(TempMax ~ fDate, data = Bats) #That is strong collinearity!
boxplot(TempMax ~ fSite, data = Bats) 
boxplot(TempMax ~ fDistRoad, data = Bats) 
boxplot(TempMax ~ fLight100, data = Bats) 
boxplot(TempMax ~ fTransect, data = Bats) 



##########################################################################


#Relationships
MyVar <- c("DayInYear", "Month", "Year",
           "DistRoad", "TempMax", 
           "WindMax")
Myxyplot(Bats, MyVar, "LACI")


xyplot(LACI ~ TempMax | Site,
       data = Bats)


plot(x = Bats$DayInYear, 
     y = Bats$LACI)
##################################################################



##################################################################
#Apply GLMM


#Check the coding of both variables
table(Bats$fDate, Bats$fTransect)
table(Bats$fSite, Bats$fTransect)
table(Bats$fDistRoad, Bats$fTransect)



#Poisson GLMM
#This is the model that was applied in the paper. 
# LACI = function of:
#          TempMax + fDistRoad + fLight100  + 
#          fSite + fYear + 
#          fSite : fDistRoad + fDistRoad : TempMax
#          random effects transect and day
#It has too many parameters for the size of the data set



#Initial analysis gave error messages for NB GLMM
#Therefore use standandized TempMax
Bats$TempMax.std <- (Bats$TempMax - mean(Bats$TempMax)) / sd(Bats$TempMax)

library(glmmADMB)
M3 <- glmmadmb(LACI ~ 1 + TempMax.std + fDistRoad + fLight100  + fSite + fYear + 
                      fSite : fDistRoad + fDistRoad : TempMax.std + 
                      (1 | fTransect) + (1 | fDate),
               data = Bats,
               family = "nbinom")
summary(M3)

######################################

#1. Set up the wind.data object

#Here are three X matrices. If you have a fast computer,
#then run the JAGS code with each of these X matrices,
#and each time write down the DIC. So..that is three JAGS
#analyses and each time use a difference X matrix. This is
#like a backwards selection, but instead of AIC we use DIC.

#Full model
X <- model.matrix(~ TempMax.std + fDistRoad + fLight100  + fSite + fYear + 
                    fSite : fDistRoad + fDistRoad : TempMax.std, 
                  data = Bats)

# Drop fSite : fDistRoad
X <- model.matrix(~ TempMax.std + fDistRoad + fLight100  + fSite + fYear + 
                    fDistRoad : TempMax.std, 
                  data = Bats)

# Drop fDistRoad : TempMax.std
 X <- model.matrix(~ TempMax.std + fDistRoad + fLight100  + fSite + fYear + 
                     fSite : fDistRoad,
                    data = Bats)

#We have already done it for you....the last one is the best
#option, as judged by the DIC.

#               				DIC
#Full model:    				803.2607
#Drop fSite : fDistRoad			802.7971
#Drop fDistRoad : TempMax.std   798.1695          
                  
#If you run the code as above, then you will be running the 
#last X matrix..which is the best model. Just go on
#running this model.
                  
K <- ncol(X) #Number of regression parameters count part

#Random effects fTransect
reTransect  <- as.numeric(as.factor(Bats$fTransect))
NumTransect <- length(unique(Bats$fTransect))

#Random effects fDate
reDate  <- as.numeric(as.factor(Bats$fDate))
NumDate <- length(unique(Bats$fDate))


#And put it all together
win.data  <- list(Y        = Bats$LACI,  #Response    
                  X        = X,          #Intercept + covariates
                  K        = K,          #Number of betas       
                  N        = nrow(Bats), #Sample size     
                  reTransect  = reTransect,  #Random effect transect
                  NumTransect = NumTransect, #Random effect transect 
                  reDate      = reDate,      #Random effect Date
                  NumDate     = NumDate      #Random effect Date 
                  )
win.data


#Step 2: Formulate JAGS modelling code
sink("BatsNBGLMM.txt")
cat("
model{
    #1A. Priors regression parameters count part and binary part
    for (i in 1:K) { beta[i]  ~ dnorm(0, 0.0001) }  
              
    #1B. Priors for random intercepts transect and date
    for (i in 1:NumTransect) {a[i] ~ dnorm(0, tau.Transect) } 
    for (i in 1: NumDate)    {b[i] ~ dnorm(0, tau.Date) } 

    #1D. Priors for variances for random intercepts
    sigma.Transect ~ dunif(0, 5)
    sigma.Date     ~ dunif(0, 5)
    tau.Transect  <- 1 / (sigma.Transect * sigma.Transect)
    tau.Date      <- 1 / (sigma.Date * sigma.Date)
    
    #1E. Half-cauchy(25) prior for size
    num.s   ~ dnorm(0, 0.0016)       #<----from NB GLM exercise
    denom.s ~ dnorm(0, 1)            #<----from NB GLM exercise
    size   <- abs(num.s / denom.s)   #<----from NB GLM exercise


    ###################################
    #2. Likelihood
    for (i in 1:N) {
       Y[i] ~  dnegbin(p[i], size)
       p[i] <- size / (size + mu[i])  
    	
       #log-link 
       log(mu[i]) <- eta[i]
       eta[i]     <- inprod(beta[], X[i,]) + 
                     a[reTransect[i]] +
                     b[reDate[i]]
                                          
       #3. Discrepancy measures: 
       #   Expected values, variance, Pearson residuals
       ExpY[i] <- mu[i] 
       VarY[i] <- mu[i] + mu[i] * mu[i] / size
       PRes[i] <- (Y[i] - ExpY[i]) / sqrt(VarY[i])   
       
       #Simulate data from a ZIP  
       YNew[i] ~  dnegbin(p[i], size) 

       #Pearson residuals for simulated data
       PResNew[i] <- (YNew[i] - ExpY[i]) / sqrt(VarY[i]) 
       D[i]       <- pow(PRes[i], 2)                      
       DNew[i]    <- pow(PResNew[i], 2)   
    }     
     Fit         <- sum(D[1:N])
     FitNew      <- sum(DNew[1:N])
}
",fill = TRUE)
sink()
#####################################
#


#Inits function
#Set initial values
inits  <- function () {
  list(beta      = rnorm(K, 0, 0.01), # betas 
       a         = rnorm(NumTransect, 0, 0.1), #random effects
       b         = rnorm(NumDate, 0, 0.1), #random effects
       sigma.Transect = runif(1, 0, 5),    #sigma transect
       sigma.Date     = runif(1, 0, 5),    #sigma Date
       num.s          = rnorm(1, 0, 25),   #for size
       denom.s        = rnorm(1, 0, 1)     #for size
       )  }


#Step 4: Parameters to estimate
params <- c("beta", "size", 
            "sigma.Transect", "sigma.Date",
            "Fit", "FitNew", "PRes",
            "ExpY"
            )

#Step 5: Run JAGS
NB1   <- jags(data        = win.data,
                inits      = inits,
                parameters = params,
                model      = "BatsNBGLMM.txt",
                n.thin     = 10, 
                n.chains   = 3,
                n.burnin   = 4000,
                n.iter     = 5000)

#And do some updating
NB2  <- update(NB1, n.iter = 10000, n.thin = 10)  
#NB3  <- update(NB2, n.iter = 50000, n.thin = 10)  
#save(NB3, file = "Bats_NB3.RData")

#It would be better to do at least another 100K iterations, 
#but this means that not all participants would be able 
#to finish the analysis.

#To reload results from an earlier session, use:
#load(file = "Bats_NB3.RData")
#recompile(NB3)
#load.module("dic")

out.NB <- NB2$BUGSoutput
# If you deciced to use NB3, then adjust the code
# above.



#Step 6: Present output 
print(out.NB, digits = 3)  


#Step 7: Assess mixing and check overdispersion
MyNames <- c(colnames(X), "sigma transect",  "sigma Day", "size")
#Adjust this variable if extra parameters are added!

MyBUGSChains(out.NB, 
             c(uNames("beta", K), "sigma.Transect", "sigma.Date", "size"),
              PanelNames = MyNames)
# We need more iterations here! But that requires half
# an hour computing time


############################################
#Check for overdispersion / underdispersion
#Bayesian approach:
mean(out.NB$sims.list$Fit >  out.NB$sims.list$FitNew) 
#That is ok-ish

#Plot the sum of squares for each iteration. 
Min <- min(out.NB$sims.list$Fit, out.NB$sims.list$FitNew)
Max <- max(out.NB$sims.list$Fit, out.NB$sims.list$FitNew)

xyplot(out.NB$sims.list$FitNew ~ out.NB$sims.list$Fit, 
       aspect = "iso", 
       col = 1, 
       pch = 16, 
       cex = 0.7,
       xlim = c(Min, Max), 
       ylim = c(Min, Max)) 
#Perfect

#What does the frequentist way of thinking tell us?
Ep <- out.NB$mean$PRes
N <- nrow(Bats)
p <- K  + 1 + 1
sum(Ep^2) / (N - p)
#Perfect


#As in the frequentist approach, we now
#need to apply a detailed model validation
#and figure out why we have overdispersion


#Step 8: Numerical output
OUT1 <- MyBUGSOutput(out.NB, 
                     c(uNames("beta", K), "sigma.Transect", "sigma.Date", "size"),
                     VarNames = MyNames)
print(OUT1, digits = 5)



#Get DIC
out.NB$DIC
#Here is where we extracted the DIC..and reran other models
#As we have already done this for you, there is no need to
#run other models.



#Task: Apply a model validation




#Calculate the fitted values.
#1. Get the betas from the MCMC iterations
#2. Define a grid of covariate values
#3. Calculate the predicted values on this 
#   grid for each MCMC iteration.
#4. Calculate the 95% credible intervals.
#5. Plot the whole thing


#1. Get the betas and gammas
beta.mcmc  <- out.NB$sims.list$beta 
dim(beta.mcmc)


#Task: make the multiplanel boxplots with model fits and 95% CIs.

