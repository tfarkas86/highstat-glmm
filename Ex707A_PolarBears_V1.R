#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  



######################################################################
#Set the working directory and import the data

setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data")
PB <- read.table(file = "PolarBearsV2.txt", 
                 header = TRUE,
                 dec = ".")
names(PB)


#Question: Is there a long term trend in the 
#          polar bear movement data?

# Model Movement as a function of Year and month
# We have multiple observations from each polar bear, hence
# use polar bear as random intercept.
# One option is:

# Movement_ij ~ N(mu_ij, sigma_Season^2) 
# E(Movement_ij) = mu_ij
# var(Movement_ij) = sigma_Season^2
# mu_ij = alpha + beta * Year + .. + Seasonality + a_i
#
#                


# But...Movement is strictly positive. We could potentially have
# negative fitted values. An alternative approach is to use a Gamma
# GLMM with log link function:
#
# Movement_ij ~ Gamma(mu_ik, r)
# E(Movement_ij) = mu_ij
# var(Movement_ij) = mu_ij^2 / r 
# mu_ij = exp(alpha + beta * Year + Seasonality + a_i)
# We can also use an identity or inverse link.



#####################################################################
#Load packages
library(lattice)
library(mgcv)
source(file = "HighstatLibV9.R")  
##################################################################


######################################################################
#Housekeeping
PB$fBearID   <- factor(PB$BearID)
PB$fRepro    <- factor(PB$Repro)
PB$fDen      <- factor(PB$Den)
PB$fMonth    <- factor(PB$Month)
PB$fSeason   <- factor(PB$Season)
######################################################################
        
        
######################################################################
#Data exploration
#Outliers
#Cleveland
plot(y = 1:nrow(PB), 
     x = PB$Movement, 
     xlab = "Values of the data",
     ylab = "Order of the data",
     pch = 16, 
     cex = 0.7)


#Relationships
plot(y = PB$Movement, 
     x = PB$dDay,
     xlab = "Days since 1 January 1988",
     ylab = "Movement")

#Add a smoother          
M1      <- loess(Movement ~ dDay, data = PB)
MyData1 <- data.frame(dDay = seq(98,4362))
P1      <- predict(M1, newdata = MyData1)
lines(x = MyData1$dDay, 
      y = P1, 
      lwd = 5)


#Month effect?
boxplot(Movement ~ fMonth,
        data = PB,
        xlab = "Month", 
        ylab = "Movement",
        varwidth = TRUE)
#Seems to be more travelling during the winter


#Movement - DayInYear effect per year
xyplot(Movement ~ DayInYear | factor(Year),
       col = 1, 
       data = PB,
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


#Movement - day effect for each bear
xyplot(Movement ~ dDay | factor(BearID),
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




#How many bears were tagged simultanously?
#Copy from here.....
MF <- function(x){ length(unique(x)) }

NBears <- tapply(PB$BearID, 
                 INDEX = PB$Year, 
                 FUN = MF)

par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 2)
plot(x = 1988:1999, 
     y = NBears, 
     pch = 16,
     xlab = "Year",
     ylab = "Number of bears with working tag")
#To here...



#Bear effect?
boxplot(Movement ~ BearID, 
        data = PB,
        xlab = "Polar bear",
        ylab = "Movement",
        varwidth = TRUE)



#Design plot
plot.design(Movement ~ fMonth + 
           fDen + fRepro + fBearID, data = PB)


#Do we have data from the same months over time?
plot(x = jitter(PB$Year), 
     y = jitter(as.numeric(PB$fMonth)),
     xlab = "Year",
     ylab = "Month",
     cex.lab = 1.5)
# Add a smoother to the graph  
M1      <- loess(as.numeric(fMonth) ~ Year, data = PB)
MyData1 <- data.frame(Year = seq(1988,1999))
P1      <- predict(M1, newdata = MyData1)
lines(x = MyData1$Year, 
      y = P1, 
      lwd = 5)

###################################################


######################################################
# Movement is strictly positive. If you apply a linear
# mixed effects model on these data, then you may end 
# up with negative fitted values and heterogeneity.
# The gamma GLM(M) avoids this.

# Let's do a little bit of theory. The code is taken
# from Chapter 6 of 'Beginner's Guide to GLM and GLMM',
# Zuur, Ieno, Hilbe (2013).


########################################################
#Section 6.5
# So..what is a GLM with a gamma distribution?

# The gamma GLM is given by:
#   Y[i]      ~ Gamma(mu[i], r)
#   E(Y[i])   = mu[i]
#   var(Y[i]) = mu[i]^2 / r
#   mu[i]     = exp(eta[i])
#   eta[i]    = Covariates

# Very similar to a negative binomial 
# GLM! The r is like the k from NB
# The only limitation is that a gamma
# is for strictly positive data: Y > 0

# Instead of a log link you can also use 
# the identity link


# So...how does the gamma GLM look like?
# COPY AND RUN THE CODE FROM HERE.....
#Figure 6.7
par(mfrow = c(2,2), mar = c(5,5,2,2))
x <- seq(0,1, length = 25)
mu <- exp(1 + 2 * x)

plot(x,mu, type = "l", cex.lab = 1.5,
     xlab = "Covariate",
     ylab = "Biomass")
x1 = seq(0,25, length = 200)

Shape <- 15
Scale <- 5 / Shape
plot(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l",
     xlab = "Possible Biomass values",
     ylab = "Probability",
     xlim = c(0,25),
     cex.lab = 1.5)

Shape <- 10
Scale <- 5 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")

Shape <- 5
Scale <- 5 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")

Shape <- 3
Scale <- 5 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")

Shape <- 2
Scale <- 5 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")

Shape <- 1
Scale <- 5 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")


###########
Shape <- 15
Scale <- 10 / Shape
plot(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l",
     xlab = "Possible Biomass values",
     ylab = "Probability",
     xlim = c(0,25),
     cex.lab = 1.5)

Shape <- 10
Scale <- 10 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")

Shape <- 5
Scale <- 10 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")

Shape <- 3
Scale <- 10 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")

Shape <- 2
Scale <- 10 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")

Shape <- 1
Scale <- 10 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")



###############
Shape <- 15
Scale <- 15 / Shape
plot(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l",
     xlab = "Possible Biomass values",
     ylab = "Probability",
     xlim = c(0,25),
     cex.lab = 1.5)

Shape <- 10
Scale <- 15 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")

Shape <- 5
Scale <- 15 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")

Shape <- 3
Scale <- 15 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")

Shape <- 2
Scale <- 15 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")

Shape <- 1
Scale <- 15 / Shape
lines(x1, dgamma(x = x1, shape = Shape, scale = Scale), type = "l")


# ALL THE WAY UP TO HERE....
# These figures show gamma distributions
# Some look like the normal distribution, others
# look completely different.




# Here is another picture what the gamma GLM
# is doing.
# COPY FROM HERE...
#Figure 6.8
par(mfrow = c(1,1))
z1 <- seq(0, 1, length = 50)
Bio <- exp(1 + 2 * z1)
Shape <- 15
Z<-matrix(nrow=50,ncol=50)
for (i in 1:50){
	Scale <- Bio[i] / Shape 
	Z[,i]<-dgamma(x=Bio, shape=Shape, scale = Scale )
}

persp(x = z1, y = Bio, z = Z,
                 scale = TRUE,
                 theta = 130, phi = 20, expand = 1,
                 ltheta = 120, shade = 0.5, 
                 ticktype = "detailed",
                 xlab = "Covariate", 
                 ylab = "Biomass", 
                 zlab = "Probability",
                 main = "")  -> res
round(res,3)
lines (trans3d(z1, y = Bio, 
               z = rep(0,50), 
               pmat = res), col = grey(0.5), lwd=5)

lines (trans3d(z1[45:50], y = Bio[45:50], 
               z = rep(0,6), 
               pmat = res), col = 1, lwd=5)

# TO HERE.....
################################################






#####################################
# Let's return to the polar bear data.
# Start of frequentist analysis


# The model that we have in mind is:
# Movement_ij = Gamma(mu_ij, r)
# E(Movement_ij) = mu_ij
# var(Movement_ij) = mu_ij ^2 / r

# log(mu_ij) = long term trend + 
#              seasonal component + 
#              covariates + polar bear effect


# But how do you do this?

# log(mu_ij) = s(Year) + s(DayInYear) + covariates + PB
# log(mu_ij) = Year + factor(Month) + covariates + PB
# log(mu_ij) = factor(Year) + factor(Season) + covariates + PB
# log(mu_ij) = Year + s(Month) + covariates + PB
# log(mu_ij) = Year + s(DayInYear) + covariates + PB
# log(mu_ij) = Year + s(DayInYear, by = factor(Year)) + covariates + PB 
# log(mu_ij) = s(DayInYear, Year) + covariates + PB 


# In this course we only discussed parametric models (no GAMMs)
# So, we can only do:

# log(mu_ij) = Year + factor(Month) + covariates + PB
#############################################


############################################
# Based on the data exploration we only use 
# fRepro and Year, and Month

# We need factors as factors
PB$fBear <- factor(PB$fBearID)
PB$fMonth <- factor(PB$Month)

#Run the Gamma glmm in glmmADMB
#It is doing the log-link.
#Computing time is about 1 minute on a fast computer.
library(glmmADMB)
M1 <- glmmadmb(Movement ~ fRepro + Year + fMonth + 
                         (1 | fBearID),
            data = PB,
            family = "gamma")
summary(M1)


# This may take too long on some computers:
# drop1(M1, test = "Chi")
              # Df   AIC   LRT  Pr(>Chi)    
# <none>           13223                    
# fRepro         2 13232 13.20   0.00136 ** 
# Year           1 13224  3.30   0.06928 .  
# factor(Month) 11 13280 79.36 1.962e-12 ***
############################################





######################
#Model validation

#Plot residuals vs each covariate
E1 <- resid(M1, type = "pearson")
E1
plot(y=1:nrow(PB), x=E1)

boxplot(E1 ~ PB$fMonth)
boxplot(E1 ~ PB$fRepro)
boxplot(E1 ~ PB$Year)

# We can calculate residuals ourselves?

# No problem...we can do that!
# Get all the components so that
# we can calculate the fitted values
# ourselves.
beta <- fixef(M1)
X    <- model.matrix(M1)
a    <- as.numeric(ranef(M1)$fBearID)
r    <- M1$alpha 
re   <- as.numeric(as.factor(PB$fBear))

# Now we are ready to calculate the fitted values,
# variance and Pearson residuals.
# Expected values of a Gamma:
mu   <- exp(X %*% beta + a[re])

# Variance of a Gamma
VarY <- mu^2 / r

# Pearson residuals:
E2 <- (PB$Movement - mu) / sqrt(VarY)
E2

# A gamma GLMM can not be overdispersed.
# And now make all the model validation graphs

# Plot residuals vs. fitted values:
plot(x = mu, y = E1)
abline(h = 0, lty = 2)


# Plot residuals vs. each covariate
plot(x = PB$Year, y = E1)
abline(h = 0, lty = 2)

boxplot(E1 ~ Month, data = PB)
abline(h = 0, lty = 2)

plot(x = PB$DayInYear , y = E1)
abline(h = 0, lty = 2)
# Check wether there is a pattern with a smoother:
library(mgcv)
T1 <- gam(E1 ~ s(DayInYear), data = PB)
summary(T1)

# There seem to be no pattern. However.....I'm worried about
# the fact that a factor(Month) term cannot capture the 
# within-month variation. That is the reason that 
# in the 'Beginner's Guide to GAM and GAMM' book we 
# used s(DayInYear).
##############################################





##############################################
#Model interpretation

#Your task:
#1. Write down the equation of the model that you just fitted.
summary(M1)
range(fitted(M1))


# Answer:
# Movement_ij     ~ Gamma(mu_ij, 1.8406)
# E(Movement_ij)   = mu_ij
# var(Movement_ij) = mu_ij^2 / 1.8406
# mu_ij = exp(Covariate stuff)


# Repro_ij = 0 & Month = 1
# mu_ij = exp(-45.30 + 0 + 0 + 0.0248 * Year_ij) 

#.....

# Repro_ij = 2 & Month = 12
# mu_ij = exp(-45.30 + 0.20 + 0.144 + 0.0248 * Year_ij) 

# Task: Why is the intercept so large (and negative)??

# Answer: Because year is from 1988! so it makes predictions for Year=0
#############################################



####################################################################
# Model interpretation
# Sketch model fit
library(plyr)
library(ggplot2)

# Create a grid with artificial covariate values.
MyData <- ddply(PB, 
                .(fRepro, fMonth), 
                summarize,
                Year = seq(min(Year), max(Year),length = 10))

X <- model.matrix(~fRepro + Year + fMonth, data = MyData)

# Get the SEs, predictor function, fitted values, and lower/upper bounds
SE  <- sqrt(diag(X %*% vcov(M1) %*% t(X)))
eta <- X %*% beta
MyData$Fit  <- exp(eta)
MyData$SeUp <- exp(eta + 2 * SE)
MyData$SeLo <- exp(eta - 2 * SE)

# And plot the whole thing
p <- ggplot()
p <- p + geom_point(data = PB, 
                    aes(y = Movement, x = Year),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Year") + 
         ylab("Movement")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = Year, 
                       y = Fit,
                       group = fMonth,
                       colour = fMonth))
#Thus gives trouble
# p <- p + geom_ribbon(data = MyData, 
                     # aes(x = Year, 
                         # ymax = SeUp, 
                         # ymin = SeLo,
                         # group = fMonth,
                         # fill = fMonth),
                     # alpha = 0.6)

p <- p + facet_grid(. ~ fRepro, scales = "fixed")
p


# Major problem:
# The SE values are way too large.
SE

# Why is this?
# Look at the first element of the following
# matrix: 
vcov(M1)


# The intercept has trouble. 
#   -What does the intercept represent? 
#   -What is the solution?
#     
#################################################





###################################################
# This is the solution:
# Standardize Year
PB$Yearc <- (PB$Year -  mean(PB$Year)) / sd(PB$Year)

#Re-apply the model
M2 <- glmmadmb(Movement ~ fRepro + Yearc + fMonth + 
                          (1 | fBearID),
              data = PB,
              family = "gamma")
summary(M2)
#Now the intercept is much smaller:

#Sketch the model again:
MyData <- ddply(PB, 
                .(fRepro, fMonth), 
                summarize,
                Yearc = seq(min(Yearc), max(Yearc),length = 10))

X   <- model.matrix(~fRepro + Yearc + fMonth, data = MyData)
SE  <- sqrt(diag(X %*% vcov(M2) %*% t(X)))
eta <- X %*% fixef(M2)
MyData$Fit  <- exp(eta)
MyData$SeUp <- exp(eta + 2 * SE)
MyData$SeLo <- exp(eta - 2 * SE)

# And plot the whole thing
# Back transform Year:
MyData$Year <- MyData$Yearc * sd(PB$Year) + mean(PB$Year)

p <- ggplot()
p <- p + geom_point(data = PB, 
                    aes(y = Movement, x = Year),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Year") + 
         ylab("Movement")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = Year, 
                       y = Fit,
                       group = fMonth,
                       colour = fMonth))

p <- p + geom_ribbon(data = MyData, 
                     aes(x = Year, 
                         ymax = SeUp, 
                         ymin = SeLo,
                         group = fMonth,
                         colour = fMonth,
                         fill = fMonth),
                     alpha = 0.1)
p <- p + facet_grid(. ~ fRepro, scales = "fixed")
p
#Not a very strong year effect!
###################################################



#######################################################
# Maybe we should have done this as the start of the 
# code?
PB$fRepro <- factor(PB$Repro)
PB$fRepro2 <- factor(PB$fRepro,
                    levels = c("0","1","2"),
                    labels = c("With cubs of the year",
                               "With one-year old cubs",
                               "With 2-year old cubs/weaning and reproducing"))
# And use Repro2 evertwhere in the code.