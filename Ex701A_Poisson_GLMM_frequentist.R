#    Introduction to Introduction to MCMC, 
#    Linear Mixed Effects models and GLMM with R
#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
#Rudy turnstone example
#Reference:
#Fuller RA, Bearhop S, Metcalfe NB, Piersma T (2013) 
#The effect of group size on vigilance in Ruddy Turnstones 
#Arenaria interpres varies with foraging habitat. 
#Ibis 155: 246-257

#4.1 Group size effect on vigilance in ruddy turnstones
#The data used in this chapter were presented in Fuller et al. (2013) 
#who investigated the effect of group size on vigilance in ruddy 
#turnstones Arenaria interpres in relation to the quality of 
#foraging habitat (in terms of food abundance and the threat of 
#their own predation) along a broadly linear 40-km stretch of 
#coastline in northeastern England.
#The ruddy turnstone is a short-legged shorebird species that 
#gets its name from its habit of turning over stones or clearing 
#away debris when it searches for intertidal food such as crustaceans, 
#worms, molluscs, and insects. 
#In their paper, Fuller et al. (2013) compared vigilance in habitat 
#types that differed greatly in prey abundance and proximity to 
#cover from which predators could launch surprise attacks.
#################################################################




#################################################################
#Set the working directory and import the data
setwd("~/Dropbox/LisbonGLMM_MCMC/")
TS <- read.table(file = "TurnstoneDataV2.txt", 
                 header = TRUE)

names(TS)
###################################################################




###################################################################
#Load packages and library files
library(lattice)  
library(lme4)
source(file = "HighstatLibV9.R")  
##################################################################


##################################################################
#Housekeeping
TS$fFlockID   <- factor(TS$FlockID)
TS$fTideState <- factor(TS$TideState)        
TS$TideStat01 <- as.numeric(TS$TideState)        
##################################################################



##################################################################
#Data exploration
#Outliers
MyVar <- c("HeadUps", "FlockSize", "NumPecks", "TimeOfDay",
           "TimeHighTide", "DaysStartWinter", "Temperature")
Mydotplot(TS[,MyVar])
#No extreme outliers

#Collinearity
MyVar <- c("FlockSize", "NumPecks", "TimeOfDay",
           "TimeHighTide", "DaysStartWinter", "Temperature",
           "TideStat01")

Mypairs(TS[,MyVar])
corvif(TS[,MyVar])

boxplot(FlockSize ~ TideState, data = TS)
boxplot(NumPecks ~ TideState, data = TS)
boxplot(TimeOfDay ~ TideState, data = TS)
boxplot(TimeHighTide ~ TideState, data = TS)
boxplot(DaysStartWinter ~ TideState, data = TS)
boxplot(Temperature ~ TideState, data = TS)

#Some trouble going on

#Relationships
MyVar <- c("FlockSize", "NumPecks", "TimeOfDay",
           "TimeHighTide", "DaysStartWinter", "Temperature")

Myxyplot(TS,MyVar,"HeadUps", MyYlab = "Head ups", MyXlab ="Covariates" )
#Some weak non-linear patterns! It is a GAMM chapter..:-)

#Balanced design?
par(mar = c(5,5,2,2))
plot(sort(table(TS$fFlockID)), 
     type = "h",
     xlab = "Flock identity",
     ylab = "Number of observations per flock",
     cex.lab = 1.5)

Z <- table(TS$fFlockID)
100 * sum(Z == 1) / length(Z)


#Zero inflation
100 * sum(TS$HeadUps == 0) / nrow(TS)



#####################################################
#Start frequentist analysis

detach(package:gamm4)
detach(package:mgcv)
detach(package:nlme)
library(lme4)

#For lme4 it may be handy to standardize covariates!
#By design of the study: Use flock id as random effect
#Based on the underlying questions and data exploration
#results we use:
M1 <- glmer(HeadUps ~ fTideState + 
                      FlockSize + 
                      NumPecks + 
                      TimeOfDay +
                      TimeHighTide  + 
                      Temperature + 
                      (1 | fFlockID),
            data = TS, 
            family = poisson)

#In the current R version, I get the following warning message:
# Warning message:
# In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
 # Model is nearly unidentifiable: very large eigenvalue
 # - Rescale variables?


#The following discussion shows a couple of things that you can try
#to avoid the warning message;

#http://stackoverflow.com/questions/23478792/warning-messages-when-trying-to-run-glmer-in-r

#The first suggestion solution is to restart from previous fitted values:
ss  <- getME(M1, c("theta","fixef"))
M1a <- update(M1, 
              start = ss, 
              control = glmerControl(optCtrl = list(maxfun=2e4)))

#Still the same warning message
#Try a different optimisation
M1b <- update(M1, 
              start = ss, 
              control = glmerControl(optimizer="bobyqa",
                               optCtrl=list(maxfun=2e4)))

# That doesn;t help neither. Try standardizing
# the covariates. Code below us not easy on the eye.

MyStd <- function(x) {  (x - mean(x))  /  sd(x) }

TS$FlockSize.c     <- MyStd(TS$FlockSize)
TS$NumPecks.c      <- MyStd(TS$NumPecks) 
TS$TimeOfDay.c     <- MyStd(TS$TimeOfDay)
TS$TimeHighTide.c  <- MyStd(TS$TimeHighTide)
TS$Temperature.c   <- MyStd(TS$Temperature)

M1c <- glmer(HeadUps ~ fTideState + 
                      FlockSize.c + 
                      NumPecks.c + 
                      TimeOfDay.c +
                      TimeHighTide.c  + 
                      Temperature.c + 
                      (1 | fFlockID),
            data = TS, 
            family = poisson)




###################################################
#Try as many different optimizers as possible
#The code below was taken from the internet. It is 
#not easy on the eye, and we will not discuss it.
#All it does is to try the GLMM with every possible
#optimizer on your computer and compare results.
#We have hashed it out, and will not discuss it. 

# afurl <- "https://raw.githubusercontent.com/lme4/lme4/master/misc/issues/allFit.R"
# library(RCurl)
# eval(parse(text=getURL(afurl)))
# aa <- allFit(M1c)
# is.OK <- sapply(aa,is,"merMod")  
# # extract the successful ones
# aa.OK <- aa[is.OK]

# #Pull out warnings:
# lapply(aa.OK,function(x) x@optinfo$conv$lme4$messages)

# #log-likelihoods are all approximately the same:
# summary(sapply(aa.OK,logLik), digits=6)

# #Compare fixed effect parameters across optimizers:
# aa.fixef <- t(sapply(aa.OK,fixef))
# library(ggplot2)
# library(reshape2)
# library(plyr)
# aa.fixef.m <- melt(aa.fixef)
# models <- levels(aa.fixef.m$Var1)
# (gplot1 <- ggplot(aa.fixef.m,aes(x=value,y=Var1,colour=Var1))+geom_point()+
    # facet_wrap(~Var2,scale="free")+
    # scale_y_discrete(breaks=models,
                     # labels=abbreviate(models,6)))


# ## coefficients of variation of fixed-effect parameter estimates:
# aa.fixef.m
# summary(unlist(daply(aa.fixef.m, "Var2", summarise, 
                     # sd(value) / abs(mean(value)))))
# #I don't think this is going good                     


# aa.varcorr <- t(sapply(aa.OK,function(x) unlist(VarCorr(x))))
# aa.varcorr.m <- melt(aa.varcorr)
# gplot1 %+% aa.varcorr.m

#End of fancy code
##########################################





#Ok..so we cannot het rid of the error...and there are 
#no smoking guns. Are the SEs ok?
summary(M1)
#No extreme large SEs...so I guess it is all ok
#Or see:
?lmerControl
#and try one of the many options


#Check overdispersion
E1 <- resid(M1, type = "pearson")
N  <- nrow(TS)
p  <- length(fixef(M1)) + 1  #'+1' is due to sigma_flock
sum(E1^2) / (N - p)
#1.18  #This is one of the very few
#examples for which a Poisson is ok!

drop1(M1c, test = "Chi")

##################################################
#Model validation M1

#Let us first calculate the fitted values.
#We can either calculate these manually, or
#use the fitted function. Let's calculate the fitted
#values manually. The model is:

# HeadUps_ij ~ Poisson(mu_ij)
# Exp(HeadUps_ij) = mu_ij
# mu_ij = exp(eta_ij)
# eta_ij = Covariate stuff + a_i
# a_i is the random intercept for flocks
# a_i ~ N(0, sigma^2_Flock)    i = 1, ..., 173

#Covariate stuff is X * beta:
X <- model.matrix( ~ fTideState + 
                     FlockSize + 
                     NumPecks + 
                     TimeOfDay +
                     TimeHighTide  + 
                     Temperature,
                     data = TS)
beta <- fixef(M1c)
CovariateStuff <- X %*% beta

#If we use X * beta as fitted values we get:
eta1 <- CovariateStuff
mu1  <- exp(eta1)   #<---- log link function!!

plot(x = mu1, 
     y =  E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline( h = 0, v = 0, lty = 2)
#looks ok

plot(x = mu1,
     y = TS$HeadUps)

#Now let's see what lme4 does
F.lme4 <- fitted(M1)

#Are these the same as our fitted values?
cbind(mu1, F.lme4) #exp(X * beta) versus exp(X * beta + random effects)
#No ....
#Why not? Because lme4 adds the random intercepts to the fitted values!!!

#This is how you extract the random effects,
#blow them up to the scale of the data,
#and calculate exp(X * beta + a)
a    <- ranef(M1)$fFlockID$'(Intercept)'

Re   <- as.numeric(as.factor(TS$fFlockID))
mu.Likelme4  <- exp(eta1 + a[Re])

cbind(mu1, mu.Likelme4, F.lme4)
#First column:  exp(X * beta)
#Second column: exp(X * beta + a)
#Third column: as second column, taken from fitted

#But if fitted() contains the random effects,
#then the resid() is also based on these values!
#Pearson residuals = (Y - fitted) / blah blah

#So..this is what lme4 will do:
E.lme4 <- (TS$HeadUps - fitted(M1) ) / sqrt(fitted(M1))
cbind(resid(M1, type = "pearson"), E.lme4)

#This is what you can do as well
E1 <- (TS$HeadUps - mu1 ) / sqrt(mu1)
cbind(E.lme4, E1)

#So..the fundamental question is whether we want to include
#the random effects in the fitted values, yes or no.
#If yes, use resid(M1) and fitted(M1)
#If no, use mu1 and E1
#But be consistent.

#At this stage it is temping to say 'yes', but the
#next exercise will show you some very odd behaviour
#if you do select 'yes'

#Let's for the moment use resid(M1) and fitted(M1)
#for model validation purposes
#So:
E.lme4 <- resid(M1, type = "pearson") #Note the type = "pearson" !!
F.lme4 <- fitted(M1)

#Plot residuals versus fitted values
plot(x = F.lme4, 
     y = E.lme4, 
     xlab = "Fitted values", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)


#Plot residuals vs each covariate
plot(x = TS$FlockSize, 
     y = E.lme4, 
     xlab = "FlockSize", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(x = TS$NumPecks, 
     y = E.lme4, 
     xlab = "NumPecks", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(x = TS$TimeOfDay, 
     y = E.lme4, 
     xlab = "TimeOfDay", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(x = TS$TimeHighTide, 
     y = E.lme4, 
     xlab = "TimeHighTide", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(x = TS$Temperature, 
     y = E.lme4, 
     xlab = "Temperature", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(x = TS$fTideState, 
     y = E.lme4, 
     xlab = "fTideState", 
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

#Step 3:
#Explain what it all means.
summary(M1c)
AIC(M1c)
drop1(M1c, test = "Chi")
#Write down the fitted equation.

# HeadUps_ij ~ Poisson(mu_ij)
# E(HeadUps_ij) = mu_ij

#High tital state:
#          0.55 - 0.016 * FlockSize + .... + 0.0063 * Temp + a_i
# mu_ij = e

#a_i ~ N(0, sigma^2_flock)


#################################################
#If you are familiar with the offset in Poisson/NB GLM(M)s,
#try using FlockSize as an offset, and explain why it is
#not a good idea to do this. 
#Offset 
TS$LFS <- log(TS$FlockSize)





#Code below is not part of this exercise
########################################
#Simulate data. 

#Let's first dump some rubbish
M2 <- glmer(HeadUps ~ fTideState + 
                      FlockSize.c + 
                      NumPecks.c + 
                      TimeOfDay.c + 
                      (1 | fFlockID),
            data = TS, 
            family = poisson)
summary(M2)


#Plan:
#Based on these results we will choose some some betas, and some
#covariate values.

#Our selected beta and sigma values.
Mybeta  <- c(1, -0.7, -0.3, 0.1, 0.25)
Mysigma <- 0.6

#Next we simulate two data sets:
#A. A nicely balanced data set
#B. A highly unbalanced data set

#Then we will extract the estimated parameters for A and B and
#compare them. And do this 1000 times


#150 clusters
#5 observations per cluster
#N = 750 observations

#Create the covariate values
I     <- sample(1:764, size = 750)
Data_750 <- TS[I,c("fTideState", "FlockSize.c", "NumPecks.c", "TimeOfDay.c")]
dim(Data_750)
head(Data_750)


#How many observations do we have?
dim(Data_750)

#Create data set A:
X <- model.matrix(~ fTideState + FlockSize.c + NumPecks.c  + TimeOfDay.c, 
                  data = Data_750)
eta   <- X %*% Mybeta
a     <- rnorm(150, mean = 0, sd = Mysigma)
Flock <- rep(1:150, each = 5)
table(factor(Flock))
mu    <- exp(eta + a[Flock])
Y     <- rpois(750, lambda = mu) 

MyData.i <- data.frame(
     HeadUps     = Y,
     fTideState  = Data_750$fTideState,
     FlockSize.c = Data_750$FlockSize.c, 
     NumPecks.c  = Data_750$NumPecks.c, 
     TimeOfDay.c = Data_750$TimeOfDay.c,
     fFlockID    = factor(Flock)
     )

M2.i <- glmer(HeadUps ~ fTideState + 
                      FlockSize.c + 
                      NumPecks.c + 
                      TimeOfDay.c + 
                      (1 | fFlockID),
            data = MyData.i, 
            family = poisson)
summary(M2.i)



#############################################
#But this is based on only one random sample
#Do it 100 times

NBoot    <- 100
Betas    <- matrix(nrow = NBoot, ncol = length(fixef(M2)))
Betas[,] <- NA
colnames(Betas) <- colnames(X)

for (i in 1: NBoot){
	print(i)
	a     <- rnorm(150, mean = 0, sd = Mysigma)
    mu    <- exp(eta + a[Flock])
    Y     <- rpois(750, lambda = mu) 

    MyData.i <- data.frame(
       HeadUps     = Y,
       fTideState  = Data_750$fTideState,
       FlockSize.c = Data_750$FlockSize.c, 
       NumPecks.c  = Data_750$NumPecks.c, 
       TimeOfDay.c = Data_750$TimeOfDay.c,
       fFlockID    = factor(Flock)
       )
  M2.i <- glmer(HeadUps ~ fTideState + 
                        FlockSize.c + 
                       NumPecks.c + 
                      TimeOfDay.c + 
                      (1 | fFlockID),
            data = MyData.i, 
            family = poisson)

  #Store regression parameters
  Betas[i,] <- fixef(M2.i)
}

Betas

#Calculate the mean and the lower/upper limites of the NBoot simulates per parameter
par2 <- colSums(Betas) / NBoot

#And the lower/upper limits
r1 <- quantile(Betas[,1], probs = c(0.025, 0.975))
r2 <- quantile(Betas[,2], probs = c(0.025, 0.975))
r3 <- quantile(Betas[,3], probs = c(0.025, 0.975))
r4 <- quantile(Betas[,4], probs = c(0.025, 0.975))
r5 <- quantile(Betas[,5], probs = c(0.025, 0.975))
SE.low <- c(r1[1], r2[1], r3[1], r4[1], r5[1])
SE.up  <- c(r1[2], r2[2], r3[2], r4[2], r5[2])


#And plot the whole things
library(coefplot2)
par1 <- fixef(M2)
se1  <- sqrt(diag(vcov(M2)))

#Observed results
coefplot2(par1, 0*se1, offset = 0, col = 1,
          varnames = names(par1), cex.var = 1,
          xlim = c(-1,1.5),
          lower2 = par1 - 1.96 * se1,
          upper2 = par1 + 1.96 * se1)

#Simulated results with 5 observations per group
coefplot2(par2, 0 * par2, offset = 0.1, col = 2,
          varnames = names(par2), 
          cex.var = 1,
          xlim = c(-1,1.5),
          lower2 = SE.low,
          upper2 = SE.up, 
          add = TRUE)



#Create data set B
#Now simulate unbalanced data (which we called B. above)

#Let's take 30% of the flocks with 1 observation
#Let's take 16% of the flocks with 2 observation
#Let's take 14% of the flocks with 3 observation
#Let's take 10% of the flocks with 4 observation
#Let's take 6% of the flocks with 5 observation
#Let's take 3% of the flocks with 6 observation
#Let's take 3% of the flocks with 7 observation
#Let's take 1% of the flocks with 8 observation

a1 <- 1:45
a2 <- rep(46:69, each = 2)
a3 <- rep(70:90, each = 3) 
a4 <- rep(91:105, each = 4)
a5 <- rep(106:114, each = 5)
a6 <- rep(115:120, each = 6)
a7 <- rep(121:127, each = 7)
a8 <- rep(128:135, each = 8)
a9 <- rep(136:144, each = 9)
a10 <- rep(145:154, each = 10)
a11 <- rep(155:165, each = 11)
a12 <- rep(166, 38)
length(c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12))


#This is our new flock id. Highly unbalanced
Flock.unb <- c(a1,a2,a3,a4,a5,a6,a7,a8,a9, a10,
              a11,a12)

table(Flock.unb)
plot(table(Flock.unb))


a     <- rnorm(166, mean = 0, sd = Mysigma)
mu    <- exp(eta + a[Flock.unb])
Y     <- rpois(750, lambda = mu) 

MyData.i <- data.frame(
     HeadUps     = Y,
     fTideState  = Data_750$fTideState,
     FlockSize.c = Data_750$FlockSize.c, 
     NumPecks.c  = Data_750$NumPecks.c, 
     TimeOfDay.c = Data_750$TimeOfDay.c,
     fFlockID    = factor(Flock.unb)
     )

M2.i <- glmer(HeadUps ~ fTideState + 
                      FlockSize.c + 
                      NumPecks.c + 
                      TimeOfDay.c + 
                      (1 | fFlockID),
            data = MyData.i, 
            family = poisson)
summary(M2.i)



#############################################
#But this is based on only one random sample
#Do it 100 times

Betas.unb    <- matrix(nrow = NBoot, ncol = length(fixef(M2)))
Betas.unb[,] <- NA
colnames(Betas.unb) <- colnames(X)

for (i in 1: NBoot){
	print(i)
	a     <- rnorm(166, mean = 0, sd = Mysigma)
    mu    <- exp(eta + a[Flock.unb])
    Y     <- rpois(750, lambda = mu) 

    MyData.i <- data.frame(
       HeadUps     = Y,
       fTideState  = Data_750$fTideState,
       FlockSize.c = Data_750$FlockSize.c, 
       NumPecks.c  = Data_750$NumPecks.c, 
       TimeOfDay.c = Data_750$TimeOfDay.c,
       fFlockID    = factor(Flock.unb)
       )
  M2.i <- glmer(HeadUps ~ fTideState + 
                        FlockSize.c + 
                       NumPecks.c + 
                      TimeOfDay.c + 
                      (1 | fFlockID),
            data = MyData.i, 
            family = poisson)

  #Store regression parameters
  Betas.unb[i,] <- fixef(M2.i)
}

Betas.unb

#Calculate the mean and the lower/upper limites of the NBoot simulates per parameter
par2.unb <- colSums(Betas.unb) / NBoot

#And the lower/upper limits
r1.unb <- quantile(Betas.unb[,1], probs = c(0.025, 0.975))
r2.unb <- quantile(Betas.unb[,2], probs = c(0.025, 0.975))
r3.unb <- quantile(Betas.unb[,3], probs = c(0.025, 0.975))
r4.unb <- quantile(Betas.unb[,4], probs = c(0.025, 0.975))
r5.unb <- quantile(Betas.unb[,5], probs = c(0.025, 0.975))
SE.low.unb <- c(r1.unb[1], r2.unb[1], r3.unb[1], r4.unb[1], r5.unb[1])
SE.up.unb  <- c(r1.unb[2], r2.unb[2], r3.unb[2], r4.unb[2], r5.unb[2])


#And plot the whole things
#Observed results
coefplot2(par1, 0*se1, offset = 0, col = 1,
          varnames = names(par1), cex.var = 1,
          xlim = c(-1,1.5),
          lower2 = par1 - 1.96 * se1,
          upper2 = par1 + 1.96 * se1)

#Simulated results with 5 observations per group
coefplot2(par2, 0 * par2, offset = 0.1, col = 2,
          varnames = names(par2), 
          cex.var = 1,
          xlim = c(-1,1.5),
          lower2 = SE.low,
          upper2 = SE.up, 
          add = TRUE)

#Simulations tudy with unbalanced dat
coefplot2(par2.unb, 0 * par2.unb, offset = 0.15, col = 3,
          varnames = names(par2.unb), 
          cex.var = 1,
          xlim = c(-1,1.5),
          lower2 = SE.low.unb,
          upper2 = SE.up.unb, 
          add = TRUE)








#Install coefplot2 (assumes shape, lme4 and coda has been installed)
# install.packages("coefplot2",
#                  repos = "http://www.math.mcmaster.ca/bolker/R",
#                  type = "source")

