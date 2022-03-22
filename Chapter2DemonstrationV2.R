#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#############################################################################


#>>> Start R code block 1 <<<

setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/")
Ospreys <- read.csv(file = "Ospreys.csv", 
                    header = TRUE, 
                    dec = ".")

str(Ospreys)
names(Ospreys)


###################################################################
#Load packages and library files
library(R2jags)
source("~/Dropbox/LisbonGLMM_MCMC//MCMCSupportHighstatV4.R")

##################################################################



#Section 2.2
M1 <- lm(THICK ~ DDD, data = Ospreys)
summary(M1)



#Figure 2.1
par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = Ospreys$DDD,
     y = Ospreys$THICK,
     xlab = "DDD",
     ylab = "Eggshell thickness")
text(0.85, 57.5, "A", cex = 1.5)

MyData <- data.frame(DDD = seq(min(Ospreys$DDD),
                               max(Ospreys$DDD),
                               length = 150))
P1 <- predict(M1, newdata = MyData, se = TRUE)

plot(x = Ospreys$DDD,
     y = Ospreys$THICK,
     xlab = "DDD",
     ylab = "Eggshell thickness")
text(0.85, 57.5, "B", cex = 1.5)


tval <- abs(qt(0.05/2, 25 - 2))
lines(MyData$DDD, P1$fit, lwd = 3)
lines(MyData$DDD, P1$fit + tval * P1$se.fit, lwd = 3, lty = 2)
lines(MyData$DDD, P1$fit - tval * P1$se.fit, lwd = 3, lty = 2)

#>>> END R code block 1 <<<
#Return to powerpoint






#>>> Start R code block 2 <<<
Ospreys$DDD.std <- scale(Ospreys$DDD)
M3 <- lm(THICK ~ DDD.std, data = Ospreys)
summary(M3)

#>>> END R code block 2 <<<
#Return to powerpoint








#>>> Start R code block 2 <<<

#Functions for the likelihood
#Will not be epxained in detail
MyLogLik <- function(theta, THICK, DDD){
	alpha <- theta[1]
	beta  <- theta[2]
	sigma <- 5.168
	N <- length(THICK)
	L <- -sum(  1 / (2*sigma*sigma)   *  
	   (THICK-alpha-beta * DDD)^2)  -  (N/2) * log(2*pi*sigma^2)
	L
}

#Function for the prior
#Will not be explained in detail
MyPrior <- function(theta){
	alpha <- theta[1]
	beta  <- theta[2]
    fprior.a <- (1/sqrt(2 * pi * 100)) * exp(-alpha^2/ (2*100))
    fprior.b <- (1/sqrt(2 * pi * 100)) * exp(-beta^2/ (2*100))
    fprior <- fprior.a * fprior.b
    log(fprior)
}

#Another support function which will not be 
#explained
JointheClub <- function(logAlpha) {
	u    <- runif(1)
    logu <- log(u)	
    logu < logAlpha
}


#Steps 1 & 2 of the MCMC algorithm
set.seed(12345)
NumberIter    <- 100000  #Get 100,000 iterations
Theta.t       <- matrix(nrow = NumberIter, ncol = 2)
Theta.star    <- vector(length = 2)
current.Theta <- c(0, 0)        #Starting values
Theta.t[1,]   <- current.Theta  #Store iteratios
#end of technical stuff

j <- 1
while (j < NumberIter){
	#EXPLAIN
	Theta.star[1] <- rnorm(1, Theta.t[j,1], 0.5)
	Theta.star[2] <- rnorm(1, Theta.t[j,2], 0.5)
	  
	#END OF R CODE BLOCK 3. RETURN TO POWERPOINT  
	#Step 3.....
	
	#Start of R code block 4
	#Run the entire algorithm
	#Calculate log(R)
	logAlpha <- MyLogLik(Theta.star, Ospreys$THICK, Ospreys$DDD.std)  + 
	            MyPrior(Theta.star) - 
	            MyLogLik(Theta.t[j,], Ospreys$THICK, Ospreys$DDD.std) - 
	            MyPrior(Theta.t[j,])
	
	
	  
	#Should the new sample join?
	Join <- JointheClub(logAlpha)
    if (Join == TRUE) { 
		j <- j + 1
		Theta.t[j,] <- Theta.star 	
    }
} 



#Theta.t contains the chains for each parameter
head(Theta.t, 10)

#Return to Powerpoint






# >>>> R code block 5 <<<
#Section 2.6.4 Mixing
par(mfrow = c(2, 1), mar = c(5, 5, 2, 2))
plot(x = 1:NumberIter,
     y = Theta.t[,1], 
     xlab = "Iterations", 
     ylab = "Intercept", 
     type = "l", 
     cex.lab = 1.5)
abline(h = coef(M3)[1], lwd =2)

plot(x = 1:NumberIter,
     y = Theta.t[,2], 
     xlab = "Iterations", 
     ylab = "Slope", 
     type = "l", 
     cex.lab = 1.5)
abline(h = coef(M3)[2], lwd =2)



#Figure 2.13
#Zoom in on the first 5000 iterations
par(mfrow = c(2, 1), mar = c(5, 5, 2, 2))
plot(x = 1:1000,
     y = Theta.t[1:1000,1], 
     xlab = "Iterations", 
     ylab = "Intercept", 
     type = "l", 
     cex.lab = 1.5)
abline(h = coef(M3)[1], lwd =2)

plot(x = 1:1000,
     y = Theta.t[1:1000,2], 
     xlab = "Iterations", 
     ylab = "Slope", 
     type = "l", 
     cex.lab = 1.5)
abline(h = coef(M3)[2], lwd =2)
#########################################





#Figure 2.14
#Chains without the burn in
BurnIn      <- 10001  #Drop the first 10,000


par(mfrow = c(2, 1), mar = c(5, 5, 2, 2))
plot(x = BurnIn:NumberIter,           #From 501 to 5000
     y = Theta.t[BurnIn:NumberIter],  #From 501 to 5000
     xlab = "Iterations", 
     ylab = "Intercept", 
     type = "l", 
     cex.lab = 1.5)
abline(h = coef(M3)[1], lwd =2)

plot(x = BurnIn:NumberIter,  
     y = Theta.t[BurnIn:NumberIter,2], 
     xlab = "Iterations", 
     ylab = "Slope", 
     type = "l", 
     cex.lab = 1.5)
abline(h = coef(M3)[2], lwd =2)

#########################################





##################################################
#Another way of plotting things
#Figure 2.15

par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = Theta.t[BurnIn:NumberIter,1],
     y = Theta.t[BurnIn:NumberIter,2],
     xlab = "Intercept",
     ylab = "Slope",
     type = "n")
abline(v = mean(Theta.t[BurnIn:NumberIter,1]), 
      col = 2, 
      lty = 2) #Posterior mean intercept
abline(h = mean(Theta.t[BurnIn:NumberIter,2]), 
       col = 2, 
       lty = 2) #Posterior mean slope

beta1 <- Theta.t[BurnIn:NumberIter,1]
beta2 <- Theta.t[BurnIn:NumberIter,2]

a0 <- c(beta1[1], beta2[1])     
for (i in 2:500) {
	 a1 <- c(beta1[i], beta2[i]) 
	 arrows(a0[1], a0[2],
	        a1[1], a1[2], length = 0)
	 a0 <- a1       
	 }
	 
	 


#Figure 2.16
###############################################################
scatterBarNorm <- function(x, dcol="blue", lhist=20, num.dnorm=5*lhist, ...){
    ## check input
    stopifnot(ncol(x)==2)
    ## set up layout and graphical parameters
    layMat <- matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
    layout(layMat, widths=c(5/7, 2/7), heights=c(2/7, 5/7))
    ospc <- 0.5 # outer space
    pext <- 4 # par extension down and to the left
    bspc <- 1 # space between scatter plot and bar plots
    par. <- par(mar=c(pext, pext, bspc, bspc),
                oma=rep(ospc, 4)) # plot parameters
    ## scatter plot
    plot(x, xlim=range(x[,1]), ylim=range(x[,2]), ...)
    ## 3) determine barplot and height parameter
    ## histogram (for barplot-ting the density)
    xhist <- hist(x[,1], plot=FALSE, breaks=seq(from=min(x[,1]), to=max(x[,1]),
                                     length.out=lhist))
    yhist <- hist(x[,2], plot=FALSE, breaks=seq(from=min(x[,2]), to=max(x[,2]),
                                     length.out=lhist)) # note: this uses probability=TRUE
    ## determine the plot range and all the things needed for the barplots and lines
    xx <- seq(min(x[,1]), max(x[,1]), length.out=num.dnorm) # evaluation points for the overlaid density
    xy <- dnorm(xx, mean=mean(x[,1]), sd=sd(x[,1])) # density points
    yx <- seq(min(x[,2]), max(x[,2]), length.out=num.dnorm)
    yy <- dnorm(yx, mean=mean(x[,2]), sd=sd(x[,2]))
    ## barplot and line for x (top)
    par(mar=c(0, pext, 0, 0))
    barplot(xhist$density, axes=FALSE, ylim=c(0, max(xhist$density, xy)),
            space=0) # barplot
    lines(seq(from=0, to=lhist-1, length.out=num.dnorm), xy, col=dcol) # line
    ## barplot and line for y (right)
    par(mar=c(pext, 0, 0, 0))
    barplot(yhist$density, axes=FALSE, xlim=c(0, max(yhist$density, yy)),
            space=0, horiz=TRUE) # barplot
    lines(yy, seq(from=0, to=lhist-1, length.out=num.dnorm), col=dcol) # line
    ## restore parameters
    par(par.)
}

X <- Theta.t[BurnIn:NumberIter,1:2]


par(mar = c(5,5,2,2))
scatterBarNorm(X,  
               xlab = "Intercept", 
               ylab = "Slope", 
               cex.lab = 1.3)







#######################################################
#Section 2.5.5 
#Present summary statistics
alpha1 <- mean(Theta.t[BurnIn:NumberIter,1], 
               na.rm = TRUE)
beta1  <- mean(Theta.t[BurnIn:NumberIter,2], 
               na.rm = TRUE)

alpha1.sd <- sd(Theta.t[BurnIn:NumberIter,1], na.rm = TRUE)
beta1.sd  <- sd(Theta.t[BurnIn:NumberIter,2], na.rm = TRUE)

#Compare home-made MCMC with glm/lm results
summary(M3)

Z <- matrix(nrow = 2, 
       c(alpha1, alpha1.sd, beta1, beta1.sd),
       byrow = TRUE)
colnames(Z) <- c("Posterior mean", "Posterior sd")
rownames(Z) <- c("Intercept", "Slope")
Z

# Small differences...probably due to the 
# small sample size

#End of R code block 5.
#Return to powerpoint



#######################################################



#R code block 6
#2.7 MCMC applied on osprey data in JAGS

#Step 1: Standardize the continuous covariates
#I am foolishly ignoring that at the moment.


#Step 2: Bundle the data
JAGS.data <- list(THICK = Ospreys$THICK, 
                  DDD   = Ospreys$DDD,
                  N     = nrow(Ospreys))
JAGS.data



#Step 3: Write JAGS code
###################################################
# 2. JAGS modelling code
sink("OspreysJAGS.txt")
cat(" 
model{
	#See powerpoint presentation on priors (4 slides)
    #1A. Priors beta and sigma
    for (i in 1:2) { beta[i] ~ dnorm(0, 0.0001)} 

    tau  <- 1 / (sigma * sigma)
    sigma ~ dunif(0, 20)


	#See powerpoint presentation on likelihood (1 slide)    
    #2. Likelihood
    for (i in 1:N) {
      THICK[i] ~ dnorm(mu[i], tau)  
      mu[i]   <- beta[1] + beta[2] * DDD[i]   
  }
 }           
",fill = TRUE)
sink()

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

#End of R block 6.
#Return to powerpoint





# >>>> Start of R code block 7  <<<<

#Step 4. Initial values 
#This is R code..not JAGS code!
inits  <- function () {
  list(
    beta  = rnorm(2, 0, 0.1),
    sigma = runif(0, 20))}


#Step 5: Parameters to save
params <- c("beta", "sigma")


#Step 6. Start JAGS
J1 <- jags(data       = JAGS.data,
           inits      = inits,
           parameters = params,
           model      = "OspreysJAGS.txt",
           n.thin     = 10,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

J2  <- update(J1, n.iter = 50000, n.thin = 10)  
out <- J2$BUGSoutput
print(out, digits = 3)  


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

#Step 7. Assess mixing
MyNames <- c("Intercept", "Slope", "sigma")

#Figure 2.18
MyBUGSChains(out, c("beta[1]", "beta[2]", "sigma"),
             PanelNames = MyNames)



#Step 8. Present output
OUT1 <- MyBUGSOutput(out, c("beta[1]"))
print(OUT1, digits = 5)

OUT1 <- MyBUGSOutput(out, c("beta[1]", "beta[2]"))
print(OUT1, digits = 5)

OUT1 <- MyBUGSOutput(out, c("beta[1]", "beta[2]", "sigma"))
print(OUT1, digits = 5)

#In case you don't want to type all the beta names:
uNames("beta", 2)
OUT1 <- MyBUGSOutput(out, c(uNames("beta", 2), "sigma"))
print(OUT1, digits = 5)


#And in case you want to modify the labels
MyNames <- c("Intercept", "Slope", "sigma")
OUT1 <- MyBUGSOutput(out, 
                     c("beta[1]", "beta[2]", "sigma"),
                     VarNames = MyNames)
print(OUT1, digits = 5)


#Posterior distribution
MyBUGSHist(out, 
           c("beta[1]", "beta[2]", "sigma"),
           PanelNames = MyNames)


#Demo code stops here. Below is some code that was
#used to make various graphs in the book chapter.
#Return to Powerpoint




#Figure 2.7
beta2 <- out$sims.matrix[,'beta[2]']
par(mar = c(5,5,2,2))
hist(beta2,
     xlab = "Posterior distribution of the slope",
     ylab = "Probability",
     main = "",
     probability = TRUE,
     cex.lab = 1.5,
     breaks = 100)

quantile(beta2, probs = c(0.25, 0.75))
quantile(beta2, probs = c(0.025, 0.975))




#Figure 2.9
beta1 <- out$sims.matrix[,'beta[1]']
beta2 <- out$sims.matrix[,'beta[2]']
sigma <- out$sims.matrix[,'sigma']
par(mfrow = c(3,1), mar = c(5,5,2,2))
hist(beta1,
     xlab = "Posterior distribution of the intercept",
     ylab = "Probability",
     main = "",
     probability = TRUE,
     cex.lab = 1.5,
     breaks = 100)

hist(beta2,
     xlab = "Posterior distribution of the slope",
     ylab = "Probability",
     main = "",
     probability = TRUE,
     cex.lab = 1.5,
     breaks = 100)

hist(sigma,
     xlab = "Posterior distribution of sigma",
     ylab = "Probability",
     main = "",
     probability = TRUE,
     cex.lab = 1.5,
     breaks = 100)




#Output used in Section 2.3.2.
m1 <- c(mean(beta1), mean(beta2), mean(sigma))
sd1 <- c(sd(beta1), sd(beta2), sd(sigma))
Z <-cbind(m1, sd1)
colnames(Z) <- c("Posterior mean", "Posterior sd")
rownames(Z) <- c("beta 1", "beta 2", "sigma")
Z
######################################################







