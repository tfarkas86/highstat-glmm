#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
#Polar bear exercise
#################################################################
#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Books/BGS/GAMM/Data/PolarBears")
PB <- read.table(file = "PolarBearsV2.txt", 
                 header = TRUE)
names(PB)

#Question: Is there a long term trend in the 
#          polar bear movement data?
# Model Movement as a function of Year and month
# We have multiple observations from each polar bear, hence
# use polar bear as random intercept.
# One option is:

# Movement_ij = alpha + beta * Year + factor(Month) + a_i + eps_ij

# This is a Gaussian mixed model with polar bear as random effect.
# But...Movement is strictly positive. We could potentially have
# negative fitted values. An alternative approach is to use a Gamma
# GLMM with log link function:
#

# Movement_ij ~ Gamma(mu_ik, r)

# E(Movement_ij) = mu_ij
# mu_ij = exp(covariate + random effects)
# var(Movement_ij) = mu_ij^2 / r 

# mu_ij = alpha + beta * Year + factor(Month) + a_i



# Different software routines have different ways of modelling
# a Gamma GLM. For example in the glm function we have:

# Movement_ij ~ Gamma(scale_ij, shape)
# E(Movement_ij)   = shape * scale_ij
# var(Movement_ij) = shape * scale_ij^2
# 
# This uses the scale and the shape. 
# And it is quite confusing to figure out what is what.
# The glm function models the scale_ij as a function of the covariates.


# This is what JAGS does:
# Y[i] ~ dgamma(r, lambda_i)
# The lamba_i in JAGS is 1 /scale_i in glm. Confusing.
# The expected value of Y[i] is r * (1/lambda_i)  

# Very confusing. In JAGS just do this:
# Y[i] ~  dgamma(r, lambda[i])
# lambda[i]  <- r / mu[i]
# log(mu[i]) <- eta[i]  
# eta[i]     <- inprod(beta[],X[i,])

#This is computer code!


#As a mathematical equation you now have:
#   Y[i]      ~ Gamma(mu[i], r)
#   E(Y[i])   = mu[i]
#   var(Y[i]) = mu[i]^2 / r
#   mu[i]     = exp(eta[i])
#   eta[i]    = X * beta + a
########################################################




########################################################
#Housekeeping
#Code factors as factors 
PB$fBearID <- factor(PB$BearID)
PB$fMonth  <- factor(PB$Month)
        
###################################################################
#Load packages and library files
library(lattice)  #Needed for multi-panel graphs
library(R2jags)
source(file="/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV9.R")  
source("/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstatV3.R")
########################################################



########################################################
#Data exploration
#Outliers
dotchart(PB$Movement, 
         xlab = "Values of the data",
         ylab = "Order of the data")




#Relationships
#Plot Movement versus time 
plot(y = PB$Movement, 
     x = PB$dDay,
     xlab = "Days since 1 January 1988",
     ylab = "Movement")
          
M1 <- loess(Movement ~ dDay, data = PB)
MyData1 <- data.frame(dDay = seq(98,4362))
P1 <- predict(M1, newdata = MyData1)
lines(MyData1$dDay, P1, lwd=5)



#Seems to be more travelling during the winter
boxplot(Movement ~ fMonth,
        data = PB,
        xlab = "Month", 
        ylab = "Movement",
        varwidth = TRUE)

xyplot(Movement ~ DayInYear | factor(Year),
       col = 1, data = PB,
       cex=0.5,
       xlab = list("Day in year", cex = 1.5),
       ylab = list("Movement", cex = 1.5),
       strip = function(bg='white',...) strip.default(bg='white', ...),
       scales = list(alternating = T,
                     x = list(relation = "same"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x,y,col=1)
       })


xyplot(Movement ~dDay | factor(BearID),
  col = 1, data = PB,
  cex=0.5,
  strip = function(bg='white',...) strip.default(bg='white', ...),
  scales = list(alternating = T,
                x = list(relation = "same"),
                y = list(relation = "same")),
  xlab = "Day number",
  ylab = "Movement",
  panel=function(x,y){
    panel.grid(h=-1, v= 2)
    panel.points(x,y,col=1)
    })


MF <- function(x){
  length(unique(x))
}


#How many bears with a working tag do we have?
NBears <- tapply(PB$BearID, INDEX = PB$Year, FUN = MF)

par(mar = c(5,5,4,2), cex.lab = 1.5)
plot(x=1988:1999, y = NBears, pch = 16, cex = 1.5,
      xlab = "Year",
      ylab = "Number of bears with working tag")


#Movement by bear
par(mar = c(5,5,4,2), cex.lab = 1.5)
boxplot(Movement ~ BearID, data = PB,
        xlab = "Polar bear",
        ylab = "Movement",
        varwidth = TRUE)


#Do we have a balanced design?
plot(x = jitter(PB$Year), 
      y = jitter(as.numeric(PB$fMonth)),
     xlab = "Year",
     ylab = "Month")
     
M1 <- loess(as.numeric(fMonth) ~ Year, data = PB)
MyData1 <- data.frame(Year = seq(1988,1999))
P1 <- predict(M1, newdata=MyData1)
lines(MyData1$Year, P1, lwd=5)
############################################################


#If you have a slow computer:
dim(PB)
Computer <- "Slow"
if (Computer == "Slow") {
	I <- sample(1:nrow(PB), size = 500, replace = FALSE)
	PB2 <- PB[I,] 
	PB2$fBearID <- factor(PB2$BearID)
    PB2$fMonth  <- factor(PB2$Month)       
    PB <- PB2  #Overwrite existing file
   }
Computer <- "Dont run again"
dim(PB)


#############################################################
#Analysis of data
#We had trouble to get a Gamma GLMM to run in glmer.

#This is a glmmADMB attempt
library(glmmADMB)
PB$Year.std      <- (PB$Year - mean(PB$Year)) / sd(PB$Year)
M1 <- glmmadmb(Movement ~ Year.std + fMonth + (1| PB$fBear), 
               family = "gamma", 
               data = PB)
summary(M1)


#And here is the JAGS code

#A full analysis of these data is given in Chapter 2
#of: A Beginner's Guide to GAMM with R (2014).
#    Zuur, Saveliev, Ieno

#We used f(DayInYear) and f(Year) as smoothers. Because we
#have not explained (and will not do so in this course) how
#to model a smoother in JAGS, we stick to parameteric models
#But a GAMM is better!



####################################

##Step 1: Fit the model.
#Use MCMC to estimate the parameters                     
#1. Bundle data
X <- model.matrix(~ Year.std + fMonth , data = PB)                      
K <- ncol(X)  
head(X)

#Random effects:
Bear <- as.numeric(as.factor(PB$fBear))
Bear
Nre <- length(unique(PB$fBear))


#The code below is copied an pasted from the owls
#Nest was changed into Plot
win.data <- list(Y    = PB$Movement, 
                 X    = X,
                 N    = nrow(PB),
                 K    = K,
                 Bear = Bear,
                 Nre  = Nre)
win.data

###################################################
# 2. JAGS modelling code
sink("GLMMPBJAGS.txt")
cat("
model{
    #1A. Priors beta and sigma
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}

    #1B. Priors random effects and sigma_Bear
    for (i in 1:Nre) { a[i] ~ dnorm(0, tau_Bear)}
   
    tau_Bear <- 1 / (sigma_Bear * sigma_Bear)
    num           ~ dnorm(0, 0.0016)    #<----half-Cauchy(25)
    denom         ~ dnorm(0, 1)         #<----half-Cauchy(25)
    sigma_Bear   <- abs(num / denom)    #<----half-Cauchy(25)

    #1C. Prior for r parameter of Gamma distribution
    r ~ dunif(0, 10)

    #2. Likelihood
    for (i in 1:N) {
      Y[i]        ~ dgamma(r, mu.eff[i])
      mu.eff[i]  <- r / mu[i]     
      log(mu[i]) <- eta[i]
      eta[i]     <- inprod(beta[], X[i,]) + a[Bear[i]] 
     }          
}
",fill = TRUE)
sink()
#####################################

#3. Initial values & parameters to save
inits  <- function () {
  list(
    beta  = rnorm(ncol(X), 0, 0.01),
    a     = rnorm(Nre, 0, 0.01),
    num   = rnorm(1, 0, 25), 
    denom = rnorm(1, 0, 1),
    r     = runif(1, 0, 10))  }

#This now becomes an issue for this data size.
#Let's take the minimal output
params <- c("beta", "a", 
            "sigma_Bear", "r")


#Don't forget to change the file name!
#4. Start JAGS
J1 <- jags(data       = win.data,
           inits      = inits,
           parameters = params,
           model      = "GLMMPBJAGS.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

#If you have a slow computer, run this:
J2  <- update(J1, n.iter = 1000, n.thin = 1)  
out <- J2$BUGSoutput

#If you have a fast computer...run this..else:
J2  <- update(J1, n.iter = 10000, n.thin = 10)  
out <- J2$BUGSoutput



print(out, digits = 3)  #Takes a few minutes
#Write down DIC at this point


#5. Assess mixing
MyBUGSChains(out, c(uNames("beta", ncol(X)), "sigma_Bear", "r"))
MyBUGSACF(out, c(uNames("beta", ncol(X)), "sigma_Bear", "r"))


#6. Present output
K <- ncol(X)
OUT1 <- MyBUGSOutput(out, c(uNames("beta", K), "sigma_Bear", "r"))
rownames(OUT1)[1:K] <- colnames(X)
print(OUT1, digits =5)
MyBUGSHist(out, c(uNames("beta", K), "sigma_Bear", "r"))
#Some parameters are not significant!
#Use DIC to figure out whether Month and Year are important



