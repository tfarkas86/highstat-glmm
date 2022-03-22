#    Highland Statistics Ltd.
#    www.highstat.com
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


####################################################


#R code Block 1

#Simulation study to understand the nature
#of Poisson and zero inflated Poisson (ZIP)
#GLMs


#Load packages
library(rgl)
library(ggplot2)
###############################################


###############################################
#First we simulate data from a Poisson GLM
#Let's assume we have 1 covariate, 100 observations,
#and the following values for the intercept and slope.
set.seed(123) #Now we all have the same results.

#Intercept, slope, sample size, covariate:
beta1 <- 1
beta2 <- 2
N     <- 1000
X1    <- runif(N)

#Predictor function eta and expected values mu
eta <- beta1 + beta2 * X1
mu  <- exp(eta)

#Observed counts:
Y   <- rpois(N, lambda = mu)

#Plot the simulated data
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = X1,
     y = Y,
     xlab = "Covariate X1",
     ylab = "Response variable Y",
     cex.lab = 1.5,
     pch = 16)

#And we can fit a Poisson GLM on these data
M1 <- glm(Y ~ X1, family = poisson)
summary(M1)


#Model validation

#Plot residuals versus fitted values
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = fitted(M1),
     y = resid(M1, type = "pearson"),
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

#End of R code block 1
#Return to Powerpoint
###############################################







##############################################
#R code Block 2
#Dispersion statistic
E1 <- resid(M1, type = "pearson")
N <- length(X1)
p <- 2
sum(E1^2) / (N - p)


#################################
#We can also sketch the model fit
#A. Create artifical covariate values
#B. Convert these into X
#C. Calculate the predictor function eta,
#   expected values ExpY, SE, and upper and
#   lower confidence intervals.
MyData      <- data.frame(X1 = seq(0, 1, length = 25))
X           <- model.matrix(~ X1, data = MyData)
MyData$eta  <- X %*% coef(M1)
MyData$ExpY <- exp(MyData$eta)
MyData$SE   <- sqrt(diag(X %*% vcov(M1) %*% t(X)))
MyData$ub   <- exp(MyData$eta + 1.96 * MyData$SE)
MyData$lb   <- exp(MyData$eta - 1.96 * MyData$SE)

#Put all relevant variables in a data frame 
PoissonData1 <- data.frame(Y  = Y,
                           X1 = X1)

p1 <- ggplot()
p1 <- p1 + geom_point(data = PoissonData1, 
                      aes(y = Y, x = X1),
                      shape = 16, 
                      size = 2)
p1 <- p1 + xlab("Covariate X1") + ylab("Response variable Y")
p1 <- p1 + theme(text = element_text(size=15))
p1 <- p1 + geom_line(data = MyData, 
                   aes(x = X1, y = ExpY), 
                   colour = I("black"),
                   lwd = 2)

p1 <- p1 + geom_ribbon(data = MyData, 
                       aes(x = X1, 
                           ymax = ub, 
                           ymin = lb ),
                       alpha = 0.2)

p1

#You may want to run the code above a couple of times
#to get a feeling about the variation in the data.
#But don't run the set.seed(12345) command again
#or else you will get the same results!
###################################################



# What is the Poisson distribution doing?
# The code below is sketching the fitted values
# again, but this time in 3-d. And the Poisson
# density curves have been added


#Three variable for 3-d scatterplot
x <- MyData$X1
y <- MyData$ExpY
z <- 0 * x

library(scatterplot3d) 
par(mfrow = c(1,1))
rr<- scatterplot3d(x, y, z, 
                   highlight.3d = FALSE, 
                   col.axis = "black",
                   col.grid = "black", 
                   pch = 20,
                   zlim = c(0, 1),
                   ylim = c(0, 30),
                   type="l",
                   lwd = 3,
                   grid = FALSE,
                   box = TRUE,
                   cex.lab = 1.5,
                   xlab = "X1",
                   ylab = "Possible values",
                   zlab = "Probability")


#Pick 5 values along the x axis to plot the density curves.
MyDatai <- data.frame(X1 = seq(from = quantile(X1, 0.15),
                                   to =   0.9 * max(X1),
                                   length = 4))


for (i in 1:4){
  Xi <- model.matrix(~ X1, data = MyDatai)[i,]
  mu.i = exp(Xi %*% coef(M1))
  yseq = round(seq(0, 30, by = 1))
  zi = dpois(yseq, lambda = mu.i)
   
  #Dotted line 
  rb = cbind(MyDatai$X1[i], yseq, 0)
  rr$points3d(rb, col = 1, type = "l", pch = ".", lty = 2)

  #points    
  rdat <- cbind(X1, Y, rep(0,N))
  rr$points3d(rdat, col = 3, type = "p", pch = 16, cex = 0.3)
 
  #Density curves 
  rb = cbind(MyDatai$X1[i], yseq, zi)
  rr$points3d(rb, col = 2, type = "h", pch = ".", lwd = 2)
 }
#




#End of R code block 2
#Return to Powerpoint
###################################################




##############################################
#R code Block 3

#Simulate Poisson data with two covariates
set.seed(123) #Now we all have the same results.

#Intercept, slope, sample size, covariate:
beta0 <- 1
beta1 <- 2
beta2 <- -3
N     <- 1000
X1    <- runif(N)
X2    <- runif(N)

#Predictor function eta and expected values mu
eta <- beta0 + beta1 * X1 + beta2 * X2
mu  <- exp(eta)

#Observed counts:
Y   <- rpois(N, lambda = mu)

#And we can fit a Poisson GLM on these data
M2 <- glm(Y ~ X1 + X2, family = poisson)
summary(M2)

#################################
#We can also sketch the model fit
#A. Create artifical covariate values
#B. Convert these into X
#C. Calculate the predictor function eta,
#   expected values ExpY, SE, and upper and
#   lower confidence intervals.

MyData      <- expand.grid(X1 = seq(0, 1, length = 25),
                           X2 = seq(0, 1, length = 25))
X           <- model.matrix(~ X1 + X2, data = MyData)
MyData$eta  <- X %*% coef(M2)
MyData$ExpY <- exp(MyData$eta)
MyData$SE   <- sqrt(diag(X %*% vcov(M2) %*% t(X)))
MyData$ub   <- exp(MyData$eta + 1.96 * MyData$SE)
MyData$lb   <- exp(MyData$eta - 1.96 * MyData$SE)

#Put all relevant variables in a data frame 
PoissonData2 <- data.frame(Y  = Y,
                           X1 = X1,
                           X2 = X2)


plot3d(x = PoissonData2$X1 ,
       y = PoissonData2$X2,
       z = PoissonData2$Y,
       type = "p",
       size = 3,
       lit = FALSE,
       xlab = "Covariate X1",
       ylab = "Covariate X2",
       zlab = "Response Y",
       col = "black")

#Add the surface for the fitted Poisson values
#For this we need to have the 25 X1 and 25 X2 values
#that we used to create the grid. Because we simulated
#them from a uniform distribution, they are nicely spread
#between 0 and 1. Later on, when we use real covariates, quite
#often multiple observations have the same covariate values. Then
#it may help the visualisation process to add some random noise.
X1.25 = seq(0, 1, length = 25)
X2.25 = seq(0, 1, length = 25)

#And we convert the vector with expected values, ExpY, into
#a 25 by 25 matrix
ExpY.2d <- matrix(MyData$ExpY, nrow = length(X1.25), ncol = length(X1.25))

#And we are ready to plot the expected values
surface3d(X1.25, X2.25, ExpY.2d, 
          alpha = 0.6, 
          front = "lines", 
          back = "lines", 
          color = "black")
#End of R code block 3
#Return to Powerpoint
###################################################





##############################################
#R code Block 4

#Bernoulli GLM
Crocs <- read.table("~/Dropbox/LisbonGLMM_MCMC/Crocodiles.txt", 
                    header = TRUE,
                    dec = ".")
str(Crocs)
names(Crocs)

#This data set is a subset of the data analyzed in Fukuda et al.
#We dropped some of the covariates.
#In the paper an information theoretic approach is presented in
#which 15-ish models are compared. Here we only focus on one model

M1 <- glm(Survived01 ~ DeltaWeight, 
          data = Crocs, 
          family = binomial(link = "logit"))
summary(M1)


#Model validation for 0-1 data is an art
E1 <- resid(M1, type = "pearson")
F1 <- fitted(M1)

par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = F1,
     y = E1,
     cex.lab = 1.5,
     ylab = "Pearson residuals",
     xlab = "Fitted probabilities")
abline(h = 0, lty = 2)     


#Model interpretation
range(Crocs$DeltaWeight)

MyData <- expand.grid(DeltaWeight = seq(from = -85.97, 
                                        to = 510.06, 
                                        length = 25))

P1          <- predict(M1, newdata = MyData, se = TRUE, type = "link")
MyData$Pi   <- exp(P1$fit) / (1 + exp(P1$fit))
MyData$SeUp <- exp(P1$fit + 1.96*P1$se.fit) / (1 + exp(P1$fit + 1.96*P1$se.fit))
MyData$SeLo <- exp(P1$fit - 1.96*P1$se.fit) / (1 + exp(P1$fit - 1.96*P1$se.fit))
MyData


p <- ggplot()
p <- p + geom_point(data = Crocs, 
                    aes(y = Survived01, x = DeltaWeight),
                    shape = 1, 
                    size = 2)
p <- p + xlab("DeltaWeight") + ylab("Probability of survival")
p <- p + theme(text = element_text(size=15))
p <- p + geom_line(data = MyData, 
                   aes(x = DeltaWeight, 
                       y = Pi), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = DeltaWeight, 
                         ymax = SeUp, 
                         ymin = SeLo ),
                     alpha = 0.2)
p





# What is the Bernoulli distribution doing?
# The code below is sketching the fitted values
# again, but this time in 3-d. And the Poisson
# density curves have been added


#Three variable for 3-d scatterplot
x <- MyData$DeltaWeight
y <- MyData$Pi
z <- 0 * x

library(scatterplot3d) 
par(mfrow = c(1,1))
rr<- scatterplot3d(x, y, z, 
                   highlight.3d = FALSE, 
                   col.axis = "black",
                   col.grid = "black", 
                   pch = 20,
                   zlim = c(0, 1),
                   ylim = c(0, 1),
                   type="l",
                   lwd = 3,
                   grid = FALSE,
                   box = TRUE,
                   cex.lab = 1.5,
                   xlab = "DeltaWeight",
                   ylab = "Possible values",
                   zlab = "Probability")


#Pick 5 values along the x axis to plot the density curves.
MyDatai <- data.frame(DeltaWeight = seq(from = quantile(Crocs$DeltaWeight, 0.15),
                                   to =   0.9 * max(Crocs$DeltaWeight),
                                   length = 4))


for (i in 1:4){
  Xi <- model.matrix(~ DeltaWeight, data = MyDatai)[i,]
  mu.i = exp(Xi %*% coef(M1)) / (1 + exp(Xi %*% coef(M1)))
  yseq = round(seq(0, 1, by = 1))
  zi = dbinom(yseq, size = 1, prob = mu.i)
   
  #Dotted line 
  rb = cbind(MyDatai$DeltaWeight[i], yseq, 0)
  rr$points3d(rb, col = 1, type = "l", pch = ".", lty = 2)

  #points    
  rdat <- cbind(Crocs$DeltaWeight, Crocs$Survived01, rep(0,nrow(Crocs)))
  rr$points3d(rdat, col = 3, type = "p", pch = 16, cex = 0.3)
 
  #Density curves 
  rb = cbind(MyDatai$DeltaWeight[i], yseq, zi)
  rr$points3d(rb, col = 2, type = "h", pch = ".", lwd = 2)
 }
#





#End of R code block 4
#Return to Powerpoint
###################################################



