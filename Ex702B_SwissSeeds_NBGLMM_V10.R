#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


################################################################
#Swiss plant seed exercise
################################################################


#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/")
SPS <- read.table("Tmontanum.txt", 
                  header = TRUE,
                  dec = ".")
##########################################################



##########################################################
# Data from
# Plant population differentiation and climate change: 
# responses of grassland species along an elevational gradient# ER FREI, J GHAZOUL, P MATTER, M HEGGLI, AR PLUESS
# Global Change Biology (2014) 20, 441–455, 
# doi: 10.1111/gcb.12403

# The aim of the study is to test for adaptation to 
# elevation of population origin, pheno- typic plasticity 
# and response to climate change of plant.

#Set up of experiment:
#1. Seeds of a specific plant species are collected from 
#   a certain number of places (= Origin) in the Swiss Alps 
#   in the summer of 2008.

SPS$Origin <- paste(SPS$Pair, SPS$Orig, sep = ".")
table(SPS$Origin)

# All observations with the same Origin value have the same
# genetic origin


#2. Seeds are put in experimental gardens (= Site)
#   at three different altitudes: 
#   600m, 1200m and 1800m above sea level
#   At each of these altitudes we have three gardens.

table(SPS$Site, SPS$Alt)

#Altitude is a treatment variable 
#Gardens (Sites) were selected at random
  
#Half of the seeds received soil treatment 1,
#the other half received soil treatment 2.

table(SPS$Site, SPS$Alt, SPS$Soil)

#We have multiple planting beds in a garden.

# To summarise:
# The experiment consisted of nine gardens in which two different 
# soil treatments (deep and shallow) were replicated in two 
# planting beds within each garden. 

# Within each planting bed, two individuals from each of 14 
# source populations of R. bulbosus, and one individual 
# from each of Orignal populations of T. montanum and 
# B. media were grown. 

table(SPS$Bed, SPS$Site, SPS$Alt, SPS$Soil)
#28 origins....28 seeds per bed.
length(table(SPS$Origin))
table(SPS$Origin, SPS$Bed)

# Immediately before transplanting the seedlings to the 
# experimental gardens, the length of the longest leaf 
# of the seedlings was measured in R. bulbosus and T. montanum, 
# and longest leaf length x number of ramets was measured in 
# B. media seedlings. 
# These proxies for the initial plant size were used as covariates 
# in the statistical analysis.


# From May to August 2010, a whole bunch of response variables 
# was measured:
#   Date of the first open flower (flowering initiation), 
#   Date of the last flower and the 
#   Date of the first mature fruit. 
#   Flowering duration 
#   Leaf and reproductive biomass 
#   Number of inflorescences (or number of flower stalks for B. media)
#       =  measure of reproduction.

#Response variables:
# BMtL10mg	vegetative growth		mg# BMS10mg  	reproductive biomass		mg# H10f	    number of inflorescences	# S10f	    number of flower stalks	# tBB10d	flowering initiation	  	DOI# tSB10d	first mature fruit	  	DOI# dB10d	    flowering duration		days

#Due to vole damage to the plants, one site 
#needs to be dropped from the analysis. This is:
#Neus

SPS2 <- SPS[SPS$Site != "Neus",]
dim(SPS)
dim(SPS2)





# To speed up the calculation time during the class, we take a random
# sample of 500 observations from the data set. 
# NOTE: the only reason for doing this is that we don't want you 
# to waste too much time waiting for your computer to finish the 
# analysis! If you have a fast computer, then set 
# IHaveSlowComputer to FALSE

IHaveSlowComputer <- TRUE
if (IHaveSlowComputer) {
  set.seed(12345)
  I1 <- sample(1:nrow(SPS2), 300)
  SPS2.smaller <- SPS2[I1, ]
  dim(SPS2.smaller)
  SPS2 <- SPS2.smaller
  dim(SPS2)
}



##############################################################
# In this exercise we will model H10f (number of inflorescences) 
# as a function of:
#   Soil treatment * altitude treatment + Initial plant size +
#   Origin effect + Site effect + Bed-within-site effect


str(SPS2)
names(SPS2)
###############################################################





###########################################################
#Load packages and library files
#Change the path of this:

library(lattice)
source("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV9.R")
source("/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstatV3.R")
######################################################################



######################################################################
#House keeping
#Create new variable with shorter/easier names
#Initial plant size
SPS2$InitSize <- SPS2$LL09a


######################################################################
#Data exploration

#Are there missing values in the relevant variables?
MyVar <- c("H10f", "InitSize", "Origin", "Site", "Bed", "Alt", "Soil")
colSums(is.na(SPS2[,MyVar]))

#Remove these rows, and only use the relevant columns
SPS3 <- SPS2[!is.na(SPS2$H10f) & !is.na(SPS2$InitSize) , MyVar]
SPS3$Site   <- factor(SPS3$Site)
SPS3$Bed    <- factor(SPS3$Bed)
SPS3$Origin <- factor(SPS3$Origin)

dim(SPS2)
dim(SPS3)

# OUTLIERS?
Myvar <- c("InitSize", "H10f")
Mydotplot(SPS3[,Myvar])
#No outliers present!


# RELATIONSHIPS?
plot(x = SPS3$InitSize,
     y = SPS3$H10f,
     xlab = "Initial plant size",
     ylab = "Number of inflorescences")

par(mfrow = c(2,3))
boxplot(H10f ~ Alt, data = SPS3)
boxplot(H10f ~ Soil, data = SPS3)
boxplot(H10f ~ Site, data = SPS3) 
boxplot(H10f ~ Origin, data = SPS3)
boxplot(H10f ~ Bed, data = SPS3)

#Zero inflation
table(SPS3$H10f)
sum(SPS3$H10f == 0) / nrow(SPS3)
#46% of the observations are equal to 0     
  

#Interactions
bwplot(H10f ~ Alt | Soil,
   strip = strip.custom(bg = 'white'),
   cex = .5, layout = c(2, 1),
   data = SPS3,
   xlab = "Altitude",
   ylab = "Number of inflorescences",
   par.settings = list(
      box.rectangle = list(col = 1),
      box.umbrella  = list(col = 1),
      plot.symbol   = list(cex = .5, col = 1)))

xyplot(H10f ~ InitSize | Alt,
       data = SPS3,
       xlab = list("Initial size", cex = 1.5),
       ylab = list("Number of inflorescences", cex = 1.5),
       layout = c(1,3),   
       strip = function(bg = 'white', ...)
       strip.default(bg = 'white', ...),
       scales = list(alternating = T,
                     x = list(relation = "same"),
                     y = list(relation = "same")),
       panel=function(x, y, subscripts){
             panel.grid(h=-1, v= 2)
             panel.loess(x, y, col = 1)
             panel.points(x, y, col = 1, pch = 16)
             })

xyplot(H10f ~ InitSize | Soil,
       data = SPS3,
       xlab = list("Initial size", cex = 1.5),
       ylab = list("Number of inflorescences", cex = 1.5),
       layout = c(1,2),   
       strip = function(bg = 'white', ...)
       strip.default(bg = 'white', ...),
       scales = list(alternating = T,
                     x = list(relation = "same"),
                     y = list(relation = "same")),
       panel=function(x, y, subscripts){
             panel.grid(h=-1, v= 2)
             panel.loess(x, y, col = 1)
             panel.points(x, y, col = 1, pch = 16)
             })


#Note the difference in observed initial size values per combination!
#Perhaps we should truncate the data?
SPS3$fAltSoil <- factor(paste(SPS3$Alt, SPS3$Soil, sep = "."))
tapply(SPS3$InitSize, FUN = range, INDEX = SPS3$fAltSoil)
tapply(SPS3$InitSize, FUN = range, INDEX = SPS3$Alt)
tapply(SPS3$InitSize, FUN = range, INDEX = SPS3$Soil)

#This is a sensible option:
SPS4 <- SPS3[SPS3$InitSize <= 2.2,]
dim(SPS3)
dim(SPS4)
############################################################



############################################################
#Start of Bayesian analysis 

# The code below carries out the following steps:
# 1. Apply Poisson GLMM
# 2. Assess overdispersion in Poisson GLMM


####################################################
#Fit the following Poisson GLMM
#
# H10f_ij    ~ P(mu_ij) 
# log(mu_ij) = InitSize + Soil * Alt + (1| Origin) + (1 | Site/Bed)




library(R2jags)


#Step 1: Collect all data for JAGS
#Standardize the continous covariates
Mystd <- function(x) {(x - mean(x)) / sd(x)}
SPS4$InitSize.c <- Mystd(SPS4$InitSize)

#Make the covariate matrix
X <- model.matrix(~1 + InitSize.c + Soil * Alt, 
                   data = SPS4)
K <- ncol(X)

#Random effects
reOrigin  <- as.numeric(as.factor(SPS4$Origin))
NumOrigin <- length(unique(SPS4$Origin))


#glmer(Y ~ X + (1 | Site/bed))  
#glmer(Y ~ X + (1 | Site) + (1 | bed) )  


reSite  <- as.numeric(as.factor(SPS4$Site))
NumSite <- length(unique(SPS4$Site))

reBed  <- as.numeric(as.factor(SPS4$Bed))
NumBed <- length(unique(SPS4$Bed))
#The random effects Bed are defined with a unique number


win.data  <- list(Y           = SPS4$H10f, #Response    
                  X           = X,         #Intercept + Soil * Alt
                  K           = K,         #Number of betas      
                  N           = nrow(SPS4),#Sample size     
                  reOrigin    = reOrigin,  #Random effect origin
                  NumOrigin   = NumOrigin, #Random effect origin 
                  reSite      = reSite,    #Random effect Site
                  NumSite     = NumSite,   #Random effect Site
                  reBed       = reBed,     #Random effect bed
                  NumBed      = NumBed)    #Random effect bed
win.data


#Step 2: Formulate JAGS modelling code
sink("PGLMM.txt")
cat("
model{
    #1A. Priors regression parameters: X * beta
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001) }  
              
    #1B. Priors for random intercepts Origin, Site and Bed
    for (i in 1:NumOrigin) {ori[i] ~ dnorm(0, tau.Or) } 
    for (i in 1:NumSite)   {sit[i] ~ dnorm(0, tau.Si) } 
    for (i in 1:NumBed)    {bed[i] ~ dnorm(0, tau.Be) } 

    #1D. Priors for variances for random intercepts ori, sit, bed
    #In case of non-mixing, change to half-Cauchy(25) or half-Cauchy(16)
    sigma.Or ~ dunif(0, 5)
    sigma.Si ~ dunif(0, 5)
    sigma.Be ~ dunif(0, 5)
    tau.Or   <- 1 / (sigma.Or * sigma.Or)
    tau.Si   <- 1 / (sigma.Si * sigma.Si)
    tau.Be   <- 1 / (sigma.Be * sigma.Be)
    
    ###################################
    #2. Likelihood
    for (i in 1:N) {
       Y[i]       ~ dpois(mu[i])
       log(mu[i]) <- eta[i]
       eta[i]     <- inprod(beta[], X[i,]) + 
                     ori[reOrigin[i]] +
                     sit[reSite[i]] +
                     bed[reBed[i]] 
                     
                      
       #3. Discrepancy measures   
        YNew[i]   ~   dpois(mu[i]) 
        expY[i]    <- mu[i] 
        varY[i]    <- mu[i]          
        PRes[i]    <- (Y[i]  - expY[i]) / sqrt(varY[i])
        PResNew[i] <- (YNew[i] - expY[i]) / sqrt(varY[i])
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
inits  <- function () {
  list(beta      = rnorm(K, 0, 0.01),
       ori       = rnorm(NumOrigin, 0, 0.1),
       sit       = rnorm(NumSite, 0, 0.1),
       bed       = rnorm(NumBed, 0, 0.1),
       sigma.Or  = runif(1, 0, 5),
       sigma.Si  = runif(1, 0, 5),
       sigma.Be  = runif(1, 0, 5)
       )  }


#Step 4: Parameters to estimate
params <- c("beta",  
            "sigma.Or", "sigma.Si", "sigma.Be",
            "Fit", "FitNew", "PRes",
            "expY", "mu")


#Step 5: Run JAGS
J1   <- jags(data       = win.data,
             inits      = inits,
             parameters = params,
             model      = "PGLMM.txt",
             n.thin     = 10, 
             n.chains   = 3,
             n.burnin   = 4000,
             n.iter     = 5000)


ptm <- proc.time() #Start the timer
J2  <- update(J1, n.iter = 10000, n.thin = 10)  
TimeTaken <- (proc.time() - ptm) / (60 * 60)  #Stop the time
save(J2, file = "J2.RData")

#It would be better to do at least another 50K iterations, but this means
#that not all participants would be able to finish the analysis.
#J3  <- update(J2, n.iter = 10000, n.thin = 10)   


#To reload results from an earlier session, use:
#load(file = "J2.RData")
#recompile(J2)
#load.module("dic")

out <- J2$BUGSoutput


#Step 6: Present output 
print(out, digits = 3)  

#Step 7: Assess mixing and check overdispersion
MyBUGSChains(out, c(uNames("beta", K), 
                           "sigma.Or", 
                           "sigma.Si", 
                           "sigma.Be"))
# We really need more iterations here! But that requires half
# an hour computing time

MyBUGSACF(out, c(uNames("beta", K)))


#Check for overdispersion / underdispersion
mean(out$sims.list$Fit >  out$sims.list$FitNew) 
#That is overdispersion!

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


#Step 8: Numerical output
OUT1 <- MyBUGSOutput(out, c(uNames("beta", K), 
                            "sigma.Or", "sigma.Si", 
                            "sigma.Be"))
rownames(OUT1)[1:K] <- colnames(X)                            
print(OUT1, digits =5)


# Write down the DIC: DIC =  1782.6
# Compare this value with the DIC of the 
# NB GAMM / ZIP GLMM / ZIP GAMM / NB GLMM
#########################################



#########################################
#Now fit a ZIP GLMM in JAGS
#Not part of this course
#########################################







#Your homework:
#Apply model validation


#So...how could we visualise the current model?
#It may help us decided what the correct model is.

#Let's calculate the posterior mean values of
#the expected values
#And then we plot these.


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

#Using these we can puzzle together the pieces.
#Note that we are not adding the random effects
#to the count part.

#2. Define a grid of covariate values
MyData <- expand.grid(InitSize.c = seq(-1.84, 2.77, length = 25),
                      Soil      = levels(SPS4$Soil),
                      Alt       = levels(SPS4$Alt))
                       
#B. Create X matrix with expand.grid
X <- model.matrix(~ InitSize.c + Soil * Alt, data = MyData)

#C. Calculate the predicted MCMC values
#   for mu
eta  <- X %*% t(beta.mcmc) 
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
p1 <- ggplot()
p1 <- p1 + geom_point(data = SPS4, 
                    aes(y = H10f, x = InitSize.c),
                    shape = 1, 
                    size = 1)
p1 <- p1 + xlab("InitSizec") + ylab("H10f")
p1 <- p1 + theme(text = element_text(size=15)) + theme_bw()
p1 <- p1 + geom_line(data = MyData2, 
                     aes(x = InitSize.c, 
                         y = mean), #This is the mu inside MyData2!
                     colour = "black")

p1 <- p1 + geom_ribbon(data = MyData2, 
                       aes(x = InitSize.c, 
                         ymax = up, 
                         ymin = lo),
                     alpha = 0.2)
p1 <- p1 + facet_grid(Soil ~ Alt, scales = "fixed")
p1

#So..these are  expected values (without the random effects) 
#with a 95% credible interval.

