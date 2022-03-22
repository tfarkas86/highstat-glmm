#    Highland Statistics Ltd.
#    www.highstat.com
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


####################################################
#Import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/FilesMCMC_GLMM/CDMCMCGLMGLMMCourse/Data")
Owls <- read.table("Owls.txt", 
                   header = TRUE)


#Check imported data
names(Owls)
str(Owls)

#This is what we have
# [1] "Nest"               "Xcoord"             "Ycoord"            
# [4] "FoodTreatment"      "SexParent"          "ArrivalTime"       
# [7] "SiblingNegotiation" "BroodSize"          "NegPerChick"       
# [10] "NCalls"          


#Task: Model SiblingNegotiation  as a function of:
#  -FoodTreatment     (factor)   
#  -SexParent         (factor)
#  -ArrivalTime       (continuous)
#  -FoodTreatment x SexParent
#  -SexParent x ArrivalTime
#  -Nest effect (but nests were selected randomly from a large group)

###################################################################
#Load packages and library files
library(lattice)  #Needed for multi-panel graphs

#source the file: HighstatLibV9.R
source(file="/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV9.R")  
#It can be dowloaded from the course website. During this course we will extend it.


##################################################################


#SiblingNegotiation is too long....use shorter name:
Owls$NCalls <- Owls$SiblingNegotiation      



##########################################
#Data exploration
boxplot(NCalls ~ Nest, data = Owls, varwidth = TRUE)
#Nest effect?


#How many nests, and how many observations per nest do we have?
length(unique(Owls$Nest))

table(Owls$Nest)
plot(table(Owls$Nest), type = "h")
# Potential trouble with 2 nests with only 4 observations


#What is the percentage of zeros?
100 * sum(Owls$NCalls == 0, na.rm = TRUE) / nrow(Owls)
plot(table(Owls$NCalls), 
     type = "h",
     cex.lab = 1.5,
     ylab = "Frequencies")


#Look at some design and interaction plots
par(mfrow = c(1,1), cex.lab = 2, mar = c(5,5,2,2))
plot.design(NCalls ~ FoodTreatment + SexParent + Nest ,
            data = Owls)

par(mfrow = c(1,1), cex.lab = 2, mar = c(5,5,2,2))            
interaction.plot(Owls$FoodTreatment, 
                 Owls$SexParent,
                 Owls$NCalls, 
                 xlab = "FoodTreatment",
                 ylab = "Sibbling negotiation")


#Plot the data as time series
library(lattice)
xyplot(NCalls ~ ArrivalTime | Nest, 
       type = "h",
       data = Owls, 
       col = 1)
  
xyplot(NCalls ~ ArrivalTime | Nest, 
       type = "h",
       subset = (FoodTreatment == "Satiated"),
       data = Owls, 
       col = 1)
quartz() #Mac command for a second graph
xyplot(NCalls ~ ArrivalTime | Nest, 
       type = "h",
       subset = (FoodTreatment == "Deprived"),
       data = Owls, 
       col = 1)



xyplot(NCalls ~ ArrivalTime | SexParent * FoodTreatment,
       data = Owls,
       ylab = "Sibbling negotiation", 
       xlab = "Arrival time (hours)",
         panel=function(x, y){
           panel.grid(h = -1, v = 2)
           panel.points(x, y, col = 1)
           panel.loess(x, y, span = 0.2, col = 1, lwd = 2)})


#Plot the spatial position of the sites
xyplot(Ycoord / 1000 ~ Xcoord / 1000, 
       aspect = "iso", 
       col = 1,
       xlab = "X coordinate",
       ylab = "Y coordinate",
       data = Owls,
       pch = 16)
#################################################



#################################################
# Frequentist analysis

# offsets fix parameter estimate to 1, so can be used to standardize by sampling
# effort where a linear 1:1 relationship is expected. Can only be used for Poisson 
# or negative binomial, because log link allows denominator of density 
# (abundance / effort) to multiplied by model. Be sure to include as a covariate
# to check that slope is very close to 1 (or just don't do this stuff at all)

# For the offset we need:
Owls$LogBroodSize <- log(Owls$BroodSize)

# Fit the following model:
# NCalls_ij ~ Poisson(mu_ij)
# mu_ij  = exp(eta_ij)
# eta_ij = Intercept + SexParent + FoodTreatment + ArrivalTime +
#          SexParent x FoodTreatment +
#          SexParent x ArrivalTime + 1 * LogBroodSize +
#          a_i
# a_i ~ N(0, sigma_Nest^2)          
# i = 1, ..., 27



# Standardize ArrivalTime to avoid warning messages in
# glmer
Owls$ArrivalTime.std <- (Owls$ArrivalTime - mean(Owls$ArrivalTime)) / 
                        sd(Owls$ArrivalTime)

library(lme4)
M1 <- glmer(NCalls ~ SexParent * FoodTreatment +
                     SexParent * ArrivalTime.std +
                     offset(LogBroodSize) + (1 | Nest),
            family = poisson,
            data = Owls)

#Check overdispersion
E1 <- resid(M1, type = "pearson")
N  <- nrow(Owls)
p  <- length(fixef(M1)) + 1
Overdisp <- sum(E1^2) / (N - p)
Overdisp

#Why do we have overdispersion?
  #A. Outliers?                  ==> Remove them?
  #B. Missing covariates or interactions?  ==> Go back or..add a latent variable 
  #C. Zero inflation?            ==> ZIP
  #D. Large variance?            ==> NB
  #E. Correlation?               ==> GLMM
  #F. Non-linear patterns        ==> GAM(M) 
  #G. Wrong link function        ==> Change it 

F1 <- fitted(M1)

#Any heterogeneity?
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)     
#Got some big residuals!


#Residual patterns?
plot(x = Owls$ArrivalTime,
     y = E1,
     xlab = "Arrival time",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)     
#Is there a non-linear pattern in here?
library(mgcv)
T1 <- gam(E1 ~ s(ArrivalTime), data = Owls)
summary(T1)
plot(T1)
abline(h = 0, lty = 2)
#Yes...there is a non-linear pattern in here.
#We should be doing a GAMM!
#But let's go on with a GLMM for the moment

MyCex <- 3 * abs(E1) / max(E1)
xyplot(Ycoord ~ Xcoord, data=Owls, pch=16, aspect="iso", 
       cex=MyCex)
# try package "inla" to incorporate spatial autocorrelation
# Bayesian incorporation of spatial correlation is very difficult

#Is eveything significant?
summary(M1)
drop1(M1, test = "Chi")
# >> But the problem is that there is overdispersion...so we cannot trust the SEs <<


#########################
#Quick and dirty solution, for overdispersed data NOT due to zero-inflation
#Observation level random intercept
Owls$Eps <- 1:nrow(Owls)
M1.eps <- glmer(NCalls ~ SexParent * FoodTreatment +
                         SexParent * ArrivalTime.std +
                         offset(LogBroodSize) + 
                         (1 | Nest) + (1 | Eps),
            family = poisson,
            data = Owls)

E1 <- resid(M1.eps, type = "pearson")
N  <- nrow(Owls)
p  <- length(fixef(M1.eps)) + 1
Overdisp <- sum(E1^2) / (N - p)
Overdisp

# Model:
# NCalls_ij ~ Poisson(mu_ij)
# E(NCalls_ij) = mu_ij
# mu_ij = exp(Intercept + Covariate stuff + a_i + eps_ij)
# a_i    ~ N(0, sigma_Nest^2)
# eps_ij ~ N(0, sigma^2)
            
            
           
summary(M1)     #Poisson GLM with random effect Nest
summary(M1.eps) #Poisson GLM with random effect Nest + olri


library(glmmADMB)
Owls$fEps <- factor(Owls$Eps)
M1.eps.admb <- glmmadmb(NCalls ~ SexParent * FoodTreatment +
                                 SexParent * ArrivalTime +
                                 + offset(LogBroodSize) + 
                                 + (1 | Nest) + (1 | fEps),
                        family = "poisson",
                        data = Owls)
summary(M1.eps.admb)

#Plot fitted versus observed values
plot(x = fitted(M1.eps.admb),
     y = Owls$NCalls,
     xlab = "Fitted values",
     ylab = "Observed data")

E1.eps <- resid(M1.eps.admb, type = "pearson")
F1.eps <- fitted(M1.eps.admb)
plot(x = F1.eps,
     y = E1.eps)



#Or do a NB GLMM (though Poisson GAMM would be better)
M2 <- glmmadmb(NCalls ~ SexParent * FoodTreatment +
                        SexParent * ArrivalTime.std +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls,
               family = "nbinom")
summary(M2)
#Check overdispersion

E2 <- resid(M2, type = "pearson")
N  <- nrow(Owls)
p  <- length(fixef(M2)) + 1 + 1  #One sigma and size/theta/k
Overdisp <- sum(E2^2) / (N - p)
Overdisp
#Ok-ish

#Model validation on the NB GLMM
F2 <- fitted(M2)

#Any heterogeneity?
plot(x = F2,
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)     
#Ok..ish? Hmmm....


#Residual patterns?
plot(x = Owls$ArrivalTime,
     y = E2,
     xlab = "Arrival time",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)     

#Is there a non-linear pattern in here?
T2 <- gam(E2 ~ s(ArrivalTime), data = Owls)
summary(T2)
plot(T2)
abline(h = 0, lty = 2)
#Yes...there is a non-linear pattern in here.
#We should be doing a GAMM!
#But let's go on with a GLMM for the moment
##############################################################



##############################################################
#Is eveything significant?
summary(M2)
drop1(M2, test = "Chi")  #Sub-models use a different k, which will throw a scale error
# Drop the SexParent: ArrivalTime.std term?

#Nested models:
# Model 1: Y = alpha + beta1 * X1 + beta2 * X2 + eps
# Model 2: Y = alpha + beta1 * X1 +     0 * X2 + eps

# M2 is nested in M1.
# H0: beta2 = 0




M3 <- glmmadmb(NCalls ~ SexParent + FoodTreatment + ArrivalTime.std +
                         SexParent : FoodTreatment + 
                         offset(LogBroodSize) + (1 | Nest),
               data = Owls,
               family = "nbinom")

drop1(M3, test = "Chi")
#Drop the interaction term


M4 <- glmmadmb(NCalls ~ SexParent + FoodTreatment + ArrivalTime.std + 
                         offset(LogBroodSize) + (1 | Nest),
               data = Owls,
               family = "nbinom")
drop1(M4, test = "Chi")
#Go on...and on....and you will end up with:



M5 <- glmmadmb(NCalls ~  FoodTreatment + ArrivalTime.std + 
                         offset(LogBroodSize) + (1 | Nest),
               data = Owls,
               family = "nbinom")
summary(M5)

####################################################################
#Task: Write down the estimated model
# Your results may be slightly different due to different
# package versions.

# Fixed effects:
                      # Estimate Std. Error t value Pr(>|z|)    
# (Intercept)            0.65154    0.09506   6.854 7.19e-12 ***
# FoodTreatmentSatiated -0.66937    0.10732  -6.237 4.45e-10 ***
# ArrivalTime.std       -0.21836    0.04816  -4.534 5.79e-06 ***


# #Deprived:
# NCalls_ij ~ NB(mu_ij, 0.86)
# E(NCalls_ij) = mu_ij
# mu_ij  = exp(eta_ij)
# eta_ij =  0.65 - 0.21 * ArrivalTime + Log(BroodSize) + a_i

# #Satiated:
# NCalls_ij ~ NB(mu_ij, 0.86)
# E(NCalls_ij) = mu_ij
# mu_ij  = exp(eta_ij)
# eta_ij = 0.65 -0.66 - 0.21 * ArrivalTime + Log(BroodSize) + a_i





####################################################################
#Sketch fitted values
#A. Specify covariate values for predictions
#B. Create X matrix with expand.grid
#C. Calculate predicted values
#D. Calculate standard errors (SE) for predicted values
#E. Plot predicted values
#F. Plot predicted values +/- 2* SE 


#A: 
range(Owls$ArrivalTime.std)
NewData <- expand.grid(ArrivalTime.std = seq(-1.59, 2.35, 
                                             length = 10) , 
                       FoodTreatment = levels(Owls$FoodTreatment))
                       
#B. Create X matrix with expand.grid
X <- model.matrix(~ FoodTreatment + ArrivalTime.std, data = NewData)


#C. Calculate predicted values
#We haven't added the offset yet!!!!
range(Owls$BroodSize)
#Calculate the fitted values for average broodsize. Don't forget
#the log!
NewData$Pred <- X %*% fixef(M5) + log(mean(Owls$BroodSize))


#D. Calculate standard errors (SE) for predicted values
#SE of fitted values are given by the square root of
#the diagonal elements of: X * cov(betas) * t(X)
NewData$SE <- sqrt(  diag(X %*%vcov(M5) %*% t(X))  )


#E. Plot predicted values
NewData

#Extract data per treatment level
Info.Sat <- NewData[NewData$FoodTreatment == "Satiated",]
Info.Dep <- NewData[NewData$FoodTreatment == "Deprived",]

plot(x = Owls$ArrivalTime.std, y = Owls$NCalls,
     col = as.numeric(Owls$FoodTreatment),
     pch = 16, 
     cex = 0.5,
     type = "n"
     )
lines(Info.Sat$ArrivalTime.std, 
      exp(Info.Sat$P), col = 2, lwd =3)
lines(Info.Sat$ArrivalTime.std, 
      exp(Info.Sat$P+2*Info.Sat$SE), col = 2, lwd = 3, lty = 2)
lines(Info.Sat$ArrivalTime.std, 
      exp(Info.Sat$P-2*Info.Sat$SE), col = 2, lwd = 3, lty = 2)

lines(Info.Dep$ArrivalTime.std, 
      exp(Info.Dep$P), col = 1, lwd = 3)
lines(Info.Dep$ArrivalTime.std, 
      exp(Info.Dep$P + 2* Info.Dep$SE), col = 1, lwd = 3, lty = 2)
lines(Info.Dep$ArrivalTime.std, 
      exp(Info.Dep$P - 2* Info.Dep$SE), col = 1, lwd = 3, lty = 2)
      
      
#############################      
#Or make a nice ggplot2 graph
library(ggplot2)

#Create artificial covariate values
library(plyr)
MyData <- ddply(Owls, .(FoodTreatment), summarize,
                ArrivalTime.std = seq(min(ArrivalTime.std), 
                                      max(ArrivalTime.std),
                                      length = 25))
MyData

#B. Create X matrix with expand.grid
X <- model.matrix(~ FoodTreatment + ArrivalTime.std, data = MyData)

#C. Calculate standard errors (SE) for predicted values
#SE of fitted values are given by the square root of
#the diagonal elements of: X * cov(betas) * t(X)

MyData$eta  <- X %*% fixef(M5)
MyData$mu  <- exp(MyData$eta + log(mean(Owls$BroodSize)))
MyData$SE   <- sqrt( diag (X %*% vcov(M5) %*% t(X) ) ) 
MyData$seup <- exp(MyData$eta + 1.96 * MyData$SE + log(mean(Owls$BroodSize)))
MyData$selo <- exp(MyData$eta - 1.96 * MyData$SE + log(mean(Owls$BroodSize)))

library(ggplot2)  
p <- ggplot()
p <- p + xlab("Arrival time") + ylab("Noises")
p <- p + theme(text = element_text(size=15))


p <- p + geom_ribbon(data = MyData,
                     aes(x = ArrivalTime.std, 
                         ymax = seup, 
                         ymin = selo),
                     fill = "red")
p <- p + geom_line(data = MyData, 
                    aes(x = ArrivalTime.std, 
                        y = MyData$mu, 
                        size = 1),    
                    col = "black")

p <- p + geom_point(data = Owls, 
                    aes(x = ArrivalTime.std, y = NCalls),
                    position = position_jitter(width = .02),
                    color = grey(0.7),
                    size = 2)
       
p <- p + facet_grid(. ~ FoodTreatment, 
                    scales = "fixed")
p <- p + theme(legend.position="none") 
p
####################################################################################




###################
#To offset or not?
M6a <- glmer(NCalls ~ FoodTreatment + ArrivalTime.std + 
                      offset(LogBroodSize) + (1 | Nest),
             family = poisson,           
             data = Owls)


M6b <- glmer(NCalls ~ FoodTreatment + ArrivalTime.std + 
                      LogBroodSize + (1 | Nest),
             family = poisson,           
             data = Owls)
summary(M6b)

M6c <- glmer(NCalls ~  FoodTreatment + ArrivalTime.std + 
                        BroodSize + (1 | Nest),
             family = poisson,          
             data = Owls)
summary(M6c)
AIC(M6a, M6b, M6c)
#But....should we not do this for a NB GLMM?
#And should we do this for the starting model?
################################################





##################################################
#In case you have glmmADMB problems...here is the same 
#analysis, but now we use glmer.nb


library(lme4)
M2 <- glmer.nb(NCalls ~ SexParent * FoodTreatment +
                        SexParent * ArrivalTime.std +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls)
summary(M2)
#Check overdispersion

E2 <- resid(M2, type = "pearson")
N  <- nrow(Owls)
p  <- length(fixef(M2)) + 1 + 1  #One sigma and size/theta/k
Overdisp <- sum(E2^2) / (N - p)
Overdisp
#Ok-ish

#Model validation on the NB GLMM
F2 <- fitted(M2)

#Any heterogeneity?
plot(x = F2,
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)     
#Ok..ish? Hmmm....


#Residual patterns?
plot(x = Owls$ArrivalTime,
     y = E2,
     xlab = "Arrival time",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)     

#Is there a non-linear pattern in here?
T2 <- gam(E2 ~ s(ArrivalTime), data = Owls)
summary(T2)
plot(T2)
abline(h = 0, lty = 2)
#Yes...there is a non-linear pattern in here.
#We should be doing a GAMM!
#But let's go on with a GLMM for the moment
##############################################################



##############################################################
#Is eveything significant?
summary(M2)
drop1(M2, test = "Chi")  #Doesn't work!!


M2a <- glmer.nb(NCalls ~ SexParent + FoodTreatment +
                        SexParent * ArrivalTime.std +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls)

M2b <- glmer.nb(NCalls ~ SexParent * FoodTreatment +
                        SexParent + ArrivalTime.std +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls)
AIC(M2, M2a, M2b)
#anova(M2, M2a)
#anova(M2, M2b)
    # df      AIC
# M2   8 3479.269
# M2a  7 3478.051
# M2b  7 3477.541  <----
#Drop SexParent x ArrivalTime.std interaction



M3 <- glmer.nb(NCalls ~ SexParent * FoodTreatment + ArrivalTime.std +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls)

M3a <- glmer.nb(NCalls ~ SexParent + FoodTreatment + ArrivalTime.std +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls)

M3b <- glmer.nb(NCalls ~ SexParent * FoodTreatment +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls)

AIC(M3, M3a, M3b)
    # df      AIC
# M3   7 3477.541
# M3a  6 3476.286  <----
# M3b  6 3495.616
#Drop  SexParent x FoodTreatment interaction


M4 <- glmer.nb(NCalls ~ SexParent + FoodTreatment + ArrivalTime.std +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls)

M4a <- glmer.nb(NCalls ~            FoodTreatment + ArrivalTime.std +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls)

M4b <- glmer.nb(NCalls ~ SexParent                + ArrivalTime.std +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls)

M4c <- glmer.nb(NCalls ~ SexParent + FoodTreatment  +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls)
AIC(M4, M4a, M4b, M4c)
   # df      AIC
# M4   6 3476.286
# M4a  5 3474.736 <----
# M4b  5 3510.353
# M4c  5 3494.201
#Drop SexParent





M5 <- glmer.nb(NCalls ~ FoodTreatment + ArrivalTime.std +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls)

M5a <- glmer.nb(NCalls ~ ArrivalTime.std +
                         offset(LogBroodSize) + (1 | Nest),
               data = Owls)

M5b <- glmer.nb(NCalls ~ FoodTreatment  +
                        offset(LogBroodSize) + (1 | Nest),
               data = Owls)
AIC(M5, M5a, M5b)
    # df      AIC
# M5   5 3474.736
# M5a  4 3509.762
# M5b  4 3492.402

#M5 is the optimal model!


#Sub-models use a different k!
summary(M5)

####################################################################
#Task: Write down the estimated model
# Your results may be slightly different due to different
# package versions.

# Fixed effects:
                      # Estimate Std. Error t value Pr(>|z|)    
# (Intercept)            0.65154    0.09506   6.854 7.19e-12 ***
# FoodTreatmentSatiated -0.66937    0.10732  -6.237 4.45e-10 ***
# ArrivalTime.std       -0.21836    0.04816  -4.534 5.79e-06 ***


# #Deprived:
# NCalls_ij ~ NB(mu_ij, 0.86)
# E(NCalls_ij) = mu_ij
# mu_ij  = exp(eta_ij)
# eta_ij =  0.65 - 0.21 * ArrivalTime + Log(BroodSize) + a_i

# #Satiated:
# NCalls_ij ~ NB(mu_ij, 0.86)
# E(NCalls_ij) = mu_ij
# mu_ij  = exp(eta_ij)
# eta_ij = 0.65 -0.66 - 0.21 * ArrivalTime + Log(BroodSize) + a_i





####################################################################
#Sketch fitted values
#A. Specify covariate values for predictions
#B. Create X matrix with expand.grid
#C. Calculate predicted values
#D. Calculate standard errors (SE) for predicted values
#E. Plot predicted values
#F. Plot predicted values +/- 2* SE 


#A: 
range(Owls$ArrivalTime.std)
NewData <- expand.grid(ArrivalTime.std = seq(-1.59, 2.35, 
                                             length = 10) , 
                       FoodTreatment = levels(Owls$FoodTreatment))
                       
#B. Create X matrix with expand.grid
X <- model.matrix(~ FoodTreatment + ArrivalTime.std, data = NewData)

#C. Calculate predicted values
#We haven't added the offset yet!!!!
range(Owls$BroodSize)
#Calculate the fitted values for average broodsize. Don't forget
#the log!
NewData$Pred <- X %*% fixef(M5) + log(mean(Owls$BroodSize))


#D. Calculate standard errors (SE) for predicted values
#SE of fitted values are given by the square root of
#the diagonal elements of: X * cov(betas) * t(X)
NewData$SE <- sqrt(  diag(X %*%vcov(M5) %*% t(X))  )


#E. Plot predicted values
NewData

#Extract data per treatment level
Info.Sat <- NewData[NewData$FoodTreatment == "Satiated",]
Info.Dep <- NewData[NewData$FoodTreatment == "Deprived",]

plot(x = Owls$ArrivalTime.std, y = Owls$NCalls,
     col = as.numeric(Owls$FoodTreatment),
     pch = 16, 
     cex = 0.5,
     type = "n"
     )
lines(Info.Sat$ArrivalTime.std, 
      exp(Info.Sat$P), col = 2, lwd =3)
lines(Info.Sat$ArrivalTime.std, 
      exp(Info.Sat$P+2*Info.Sat$SE), col = 2, lwd = 3, lty = 2)
lines(Info.Sat$ArrivalTime.std, 
      exp(Info.Sat$P-2*Info.Sat$SE), col = 2, lwd = 3, lty = 2)

lines(Info.Dep$ArrivalTime.std, 
      exp(Info.Dep$P), col = 1, lwd = 3)
lines(Info.Dep$ArrivalTime.std, 
      exp(Info.Dep$P + 2* Info.Dep$SE), col = 1, lwd = 3, lty = 2)
lines(Info.Dep$ArrivalTime.std, 
      exp(Info.Dep$P - 2* Info.Dep$SE), col = 1, lwd = 3, lty = 2)
      
      
#############################      
#Or make a nice ggplot2 graph
library(ggplot2)

#Create artificial covariate values
library(plyr)
MyData <- ddply(Owls, .(FoodTreatment), summarize,
                ArrivalTime.std = seq(min(ArrivalTime.std), 
                                      max(ArrivalTime.std),
                                      length = 25))
MyData

#B. Create X matrix with expand.grid
X <- model.matrix(~ FoodTreatment + ArrivalTime.std, data = MyData)

#C. Calculate standard errors (SE) for predicted values
#SE of fitted values are given by the square root of
#the diagonal elements of: X * cov(betas) * t(X)

MyData$eta  <- X %*% fixef(M5)
MyData$mu  <- exp(MyData$eta + log(mean(Owls$BroodSize)))
MyData$SE   <- sqrt( diag (X %*% vcov(M5) %*% t(X) ) ) 
MyData$seup <- exp(MyData$eta + 1.96 * MyData$SE + log(mean(Owls$BroodSize)))
MyData$selo <- exp(MyData$eta - 1.96 * MyData$SE + log(mean(Owls$BroodSize)))

library(ggplot2)  
p <- ggplot()
p <- p + xlab("Arrival time") + ylab("Noises")
p <- p + theme(text = element_text(size=15))


p <- p + geom_ribbon(data = MyData,
                     aes(x = ArrivalTime.std, 
                         ymax = seup, 
                         ymin = selo),
                     fill = "red")
p <- p + geom_line(data = MyData, 
                    aes(x = ArrivalTime.std, 
                        y = MyData$mu, 
                        size = 1),    
                    col = "black")

p <- p + geom_point(data = Owls, 
                    aes(x = ArrivalTime.std, y = NCalls),
                    position = position_jitter(width = .02),
                    color = grey(0.7),
                    size = 2)
       
p <- p + facet_grid(. ~ FoodTreatment, 
                    scales = "fixed")
p <- p + theme(legend.position="none") 
p


###########################################
#3D graph

library(rgl)


#And use red for points with a larger Pearson residual
plot3d(x = Owls$ArrivalTime.std,
       y = Owls$LogBroodSize,
       z = Owls$NCalls,
       type = "p",
       size = 2,
       lit = FALSE,
       xlab = "ArrivalTime.std",
       ylab = "LogBroodSize",
       zlab = "NCalls")

#Now we would like to add the fitted values
#of the NB GLMM as a surface.
MyData <- ddply(Owls, .(FoodTreatment), summarize,
                ArrivalTime.std = seq(min(ArrivalTime.std), 
                                      max(ArrivalTime.std),
                                      length = 25))
MyData


#Calculate the fitted values
range(Owls$ArrivalTime.std)
range(Owls$LogBroodSize)
MyData <- expand.grid(ArrivalTime.std = seq(-1.59, 2.35, length = 25),
                      LogBroodSize = seq(0, 1.95, length = 25),
                      FoodTreatment = levels(Owls$FoodTreatment))
 
X    <- model.matrix(~ 1 + FoodTreatment + ArrivalTime.std, data = MyData)
eta  <- X %*% fixef(M5) + MyData$LogBroodSize
mu   <- exp(eta)
ExpY <- mu  #Expected values

#And add the expected values to MyData
MyData2 <- cbind(MyData, ExpY)

#So..these are the two variables that we use
#to calculate our grid for:
X1 <- seq(-1.59, 2.35, length = 25)
X2 <- seq(0, 1.95, length = 25)

#And we convert the vector with expected values into
#a 25 by 25 matrix
#But...we have two planes. Extract the fitted values per FoodTreatment level
ExpY.Deprived <- MyData2$ExpY[MyData2$FoodTreatment == "Deprived"]
ExpY.2d.Deprived <- matrix(ExpY.Deprived, nrow = 25, ncol = 25)

ExpY.Satiated <- MyData2$ExpY[MyData2$FoodTreatment == "Satiated"]
ExpY.2d.Satiated <- matrix(ExpY.Satiated, nrow = 25, ncol = 25)


#Add the surface for the fitted NB model
surface3d(X1, X2, ExpY.2d.Deprived, 
          alpha = 0.6, 
          front = "lines", 
          back = "lines", 
          color = "black")

surface3d(X1, X2, ExpY.2d.Satiated, 
          alpha = 0.6, 
          front = "lines", 
          back = "lines", 
          color = "red")
          
          