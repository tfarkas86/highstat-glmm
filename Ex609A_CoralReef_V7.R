
#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
#Coral reef exercise
#These data are taken from:
#Caribbean-wide decline in carbonate production threatens 
#coral reef growth
#Perry et al. (2012)
#NATURE COMMUNICATIONS | 4:1402 | DOI: 10.1038/ncomms2409 | www.nature.com/naturecommunications


#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Books/BGS/GAMM/Data/ReefData")
CO <- read.table("CoralData.txt", header = TRUE)

names(CO)
str(CO)
# 'data.frame':	101 obs. of  10 variables:
 # $ Country         : Factor w/ 4 levels "Bahamas","Belize",..: 1 1 1 1 1  ...
 # $ Site            : Factor w/ 14 levels "BAR","CAI","CAL",..: 5 5 5 7 7 7 1 ...
 # $ Transect        : int  1 2 3 1 2 3 1 2 3 1 ...
 # $ Depth           : int  19 19 19 17 17 17 20 20 20 5 ...
 # $ Habitat         : Factor w/ 5 levels "APZ","FRMZ","FRS",..: 4 4 4 4 4  ...
 # $ LCC             : num  7.88 9.61 5.69 5.11 4.09 ...
 # $ Gross_production: num  2.04 1.94 1.25 1.1 1.06 1.33 1.8 1.37 2.56 2.09 ...
 # $ Gross_erosion   : num  1.49 1.51 1.54 1.7 1.71 1.64 1.2 1.17 1.26 1.93 ...
 # $ Net_production  : num  0.56 0.43 -0.29 -0.6 -0.65 -0.31 0.6 0.2 1.3 0.16 ...
 # $ Accretion_rate  : num  0.57 0.49 0.05 -0.21 -0.23 -0.01 0.57 0.29 1.05 0.42 ... 


###################################################################
#Load packages and library files
library(lattice)  #Needed for multi-panel graphs
library(nlme)
library(mgcv)

#source the file: HighstatLibV7.R
source("~/Dropbox/LisbonGLMM_MCMC/HighstatLibV9.R")
##################################################################

#The variables:
# [1] "Country"          "Site"             "Transect"         "Depth"
# [5] "Habitat"          "LCC"              "Gross_production" "Gross_erosion"
# [9] "Net_production"   "Accretion_rate"

# Aim of the study:
# Perry et al. (2013) sampled reef carbonate 
# production (CaCO3 m–2 yr–1) at 101 transects 
# in 19 coral reefs occurring in four Caribbean 
# countries (Bahamas, Belize, Bonaire, and Grand Cayman). 
# Surveys were conducted within a range of common 
# Caribbean reef habitats: nearshore hardgrounds, 
# Acropora palmata habitats, Montastraea spur-and-groove zones, 
# fore-reef slopes, and deeper (18–20 m) shelf-edge Montastraea reefs. 
# See Perry et al. (2013) for further details.


# The response variable is Net_production, 
# and we have one observation within a transect. 
# Hence, transect is the statistical sampling unit. 
# Potential explanatory variables are Country, 
# Site (a group of transects), Depth, Habitat, and 
# LCC (life coral cover %). 

# Task:
# Model net carbonate production as a function of:
# LCC ( Live coral cover)
# Habitat
# Random effect site nested within a Country
###################################################




###################################################
#Housekeeping
CO$fCountry <- factor(CO$Country)
CO$fSite    <- factor(CO$Site)
CO$fHabitat <- factor(CO$Habitat)

#I'm lazy....this is the response variable
CO$G        <- CO$Net_production
###################################################



###################################################
#Data exploration
#Outliers
MyVar <- c("G", "LCC")
Mydotplot(CO[,MyVar])


#Are the data balanced?
table(CO$fCountry)
table(CO$fSite)
table(CO$fHabitat)

#4 countries
#14 uniquely defined sites


#Relationships
boxplot(G ~ fSite, data = CO, xlab = "Site")
boxplot(G ~ fCountry, data = CO, xlab = "Country")
boxplot(G ~ fHabitat, data = CO, xlab = "Habitat")

MyYlab <- expression(paste("Net carbonate production ", 
                "(kg CaCO"[3],
                " m"^"-2",
                "year" ^"-1",")"))
par(mar = c(5,6,2,2), cex.lab = 1.5)
plot(x = CO$LCC, 
     y = CO$G,
     pch = 16,
     cex = 1,
     xlab = "LCC (%)",
     ylab = MyYlab)
abline(h = 0, lty = 2)


#############
#Collinearity
boxplot(LCC ~ fSite, data = CO, xlab = "Site", ylab = "LCC (%)")
boxplot(LCC ~ fCountry, data = CO, xlab = "Country", ylab = "LCC (%)")

par(cex.lab = 1.5, mar = c(5,5,2,2))
boxplot(LCC ~ fHabitat, data = CO, xlab = "Habitat", ylab = "LCC (%)")

#Some strong collinearity going on here!

#There is more trouble:
table(CO$fCountry, CO$fHabitat)
#Use them both an you may get trouble!
#Some levels of country cannot be discriminated from some habitat levels.
#Options: 
# 1. Remove Bahamas and SEMR 
# 2. Drop one covariate and be careful with the interpretation of the results
# 3. Go back and get more data 

# Some data phishing indicates that country is not significant...but 
# habitat level is!
# Try presenting that in a paper!


#Data exploration indicates that there are no obvious 
#outliers; there is a positive, potentially nonlinear 
#LCC effect; there is heterogeneity; and there is a 
#certain degree of collinearity between LCC and the 
#categorical covariates. There is also trouble with country and habitat.
##########################################################





#Start of analysis. First the frequentist analysis
##########################################################
#Frequentist analysis

#Based on the Perry et al. (2012) paper, and based
#on the data exploration we start as follows:
M0 <- lme(G ~ LCC + fHabitat,
          random =~ 1 |fSite,
          data = CO, 
          method = "REML")
summary(M0)

Anova(M0, type=2)
resvar <- 1.575749; randvar <- 1.084636
icc <- randvar^2/(resvar^2 + randvar^2) # 0.321

#By design of the study we use site as random effect.
#There are only 4 countries, so we can't use a 2-way nested structure.
                   
                       
#Model validation                       
par(mfrow = c(2,2), mar = c(5,5,2,2), cex.lab = 1.4)
E0 <- resid(M0, type = "n")
plot(x = CO$LCC, 
     y = E0, 
     xlab = "LCC",
     ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(E0 ~ fCountry,   
        ylab = "Normalized residuals",
        data = CO, xlab = "Country")
abline(h = 0, lty = 2)

boxplot(E0 ~ fHabitat, 
        ylab = "Normalized residuals", 
        data = CO,  xlab = "Habitat")
abline(h = 0, lty = 2)

boxplot(E0 ~ fSite, 
         ylab = "Normalized residuals",
         data = CO, 
         xlab = "Site",
         varwidth = TRUE)
abline(h = 0, lty = 2)

#Ouch....do we see a non-linear effect of LCC? I don't think so.
#And heterogeneity? Big time.

#If a plot of the residuals versus a covariate shows 
#any patterns then the model needs to be improved. 
#Instead of plotting the residuals versus a covariate 
#it is also an option to apply a regression model 
#or GAM in which residuals are modelled as a function 
#of the covariate. A significant effect means that the 
#original model is not good. 
#Assuming you are familair with GAM, we use:

L1 <- gam(E0 ~ s(LCC), data = CO)
summary(L1)

#17% of the variation in the residuals (!) is explained
#by a smoother of LCC. And this is how the smoother looks like:

par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(L1, xlab = "LCC", cex.lab = 1.5, ylim = c(-2.5, 4))
points(x = CO$LCC, 
       y = E0)
abline(h = 0)
#Ok...there are actually only 4 points that cause the non-linear effect!

#If you are not familair with GAM, then consider it as a moving average smoother


#In Zuur, Saveliev, Ieno (2014). A Beginner's Guide to GAMM
#we use these data to show how additive models work.


#There is also hetereogeneity. There are two options to solve this:
#1. The good way: Add a variance structure to model the heterogeneity
#   See Chapter 4 in Zuur et al. (2009)..or Zuur et al. (2014)

#2. The bad way: A transformation on G. Why bad? A transformation may
#   change relationships. See our introductory course.

#But....a log transformation does make life much easier as it may
#also linearize the relationships. And Pery et al. (2012) did the
#same (not that this justifies a transformation)

#So...let's do things quick and dirty and transform the response variable
#The alternative is to come back and do a GAMM course!


# Log transformation: log(Y + constant)   
# where the constant is small and Y + constant > 0
min(CO$G)
CO$LogG <- log(CO$G + 3)  #I really prefer a GAMM here! 


#Refit the model on log transformed values
M1 <- lme(LogG ~ LCC + fHabitat,
          random =~ 1 |fSite,
          data = CO, 
          method = "REML")

#Model validation                       
par(mfrow = c(2,2), mar = c(5,5,2,2), cex.lab = 1.4)
E1 <- resid(M1, type = "n")
plot(x = CO$LCC, 
     y = E1, 
     xlab = "LCC",
     ylab = "Normalized residuals")
abline(h = 0, lty = 2)

boxplot(E1 ~ fCountry,   
        ylab = "Normalized residuals",
        data = CO, xlab = "Country")
abline(h = 0, lty = 2)

boxplot(E1 ~ fHabitat, 
        ylab = "Normalized residuals", 
        data = CO,  xlab = "Habitat")
abline(h = 0, lty = 2)

boxplot(E1 ~ fSite, 
         ylab = "Normalized residuals",
         data = CO, xlab = "Site")
abline(h = 0, lty = 2)

L1 <- gam(E1 ~ s(LCC), data = CO)
summary(L1)
plot(L1)

#Heterogeneity mostly gone...non-linear pattern with LCC certainly gone!
#But....we have essentially taken a hammer and removed some
#interesting information!




#############################################################




#Step 2: Are all terms significant?
M2 <- lme(LogG ~ LCC + fHabitat,
          random =~ 1 |fSite,
          data = CO, 
          method = "ML")
M2A <- update(M2, .~. - LCC)
M2B <- update(M2, .~. - fHabitat)
anova(M2, M2A)
anova(M2, M2B)  

#Everything significant (though we don't know yet whether
#the habitat effect is in fact a country effect)

#Step 3: Explain what it all means.
summary(M1)  #This is REML obtained stuff
#Calculate the intra-class correlation.
0.1563258 ^2 / (0.1563258 ^2 +  0.2795533 ^2) # 0.238

#Write down the 5 equations for LogG (one per habitat type).
# LogG_ij  ~ N(mu_ij, sigma^2)
# E(LogG_ij) = mu_ij
# mu_ij      = Covariate effects + a_i + eps_ij
# a_i    ~ N(0, sigma_Site^2)
# eps_ij ~ N(0, sigma^2)
# where i = 1, ..., 14

#The expression for mu_ij is as follows.
#APZ:
#mu_ij = 0.432 + 0.048 * LCC_ij + a_i 

#FRMZ
#mu_ij = 0.432 + 0.369 + 0.048 * LCC_ij + a_i 
#      = 0.807 + 0.048 * LCC_ij + a_i 
#....

#SHG:
#mu_ij = 0.432 - 0.246 + 0.048 * LCC_ij + a_i 
#      = 0.192 + 0.048 * LCC_ij + a_i



##############################################
#Sketch fitted values
#A. Specify covariate values for predictions
#B. Create X matrix with expand.grid
#C. Calculate predicted values
#D. Calculate standard errors (SE) for predicted values
#E. Plot predicted values
#F. Plot predicted values +/- 	1.96 * SE


MyData <- expand.grid(LCC = seq(0.67, 48.35, length = 10),
                      fHabitat = levels(CO$fHabitat))

#Or better:
#Here is some fancy code to make a better MyData
library(plyr)
MyData <- ddply(CO, .(fHabitat), summarize,
                LCC = seq(min(LCC), max(LCC),length = 10))
MyData



X  <- model.matrix(~ LCC + fHabitat, data = MyData)
mu <- X %*% fixef(M2)
SE <- sqrt(diag(X %*% vcov(M2) %*% t(X)))


#Glue them all together
MyData$mu <- mu
MyData$ub <- mu + 1.96 * SE
MyData$lb <- mu - 1.96 * SE
head(MyData)


library(ggplot2)
p <- ggplot()
p <- p + geom_point(data = CO, 
                    aes(y = LogG, x = LCC),
                    shape = 16, 
                    size = 3)
p <- p + xlab("LCC") + ylab("LogG")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = LCC, y = mu), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = LCC, 
                         ymax = ub, 
                         ymin = lb ),
                     alpha = 0.5)
p <- p + facet_grid(. ~ fHabitat, scales = "fixed")
p
#How can you criticise this graph?


##########################################################
#And tidy up this graph!
tapply(CO$LCC, FUN = range, INDEX = CO$fHabitat)                     

LCC1 <- seq(2.37, 48.35, length = 25)         
LCC2 <- seq(4.06, 45.75, length = 25)         
LCC3 <- seq(3.91, 31.02, length = 25)         
LCC4 <- seq(4.09, 12.51, length = 25)         
LCC5 <- seq(0.67, 10.16, length = 25)         

MyData <- data.frame(LCC = c(LCC1, LCC2, LCC3, LCC4, LCC5),
                     fHabitat = rep(levels(CO$fHabitat), each = 25))
                     
X      <- model.matrix(~ LCC + fHabitat, data = MyData)

mu <- X %*% fixef(M2)
SE <- sqrt(diag(X %*% vcov(M2) %*% t(X)))


#Glue them all together
MyData$mu <- mu
MyData$ub <- mu + 1.96 * SE
MyData$lb <- mu - 1.96 * SE
MyData

p <- ggplot()
p <- p + geom_point(data = CO, 
                    aes(y = LogG, x = LCC),
                    shape = 16, 
                    size = 3)
p <- p + xlab("LCC") + ylab("LogG")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = LCC, y = mu), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = LCC, 
                         ymax = ub, 
                         ymin = lb ),
                     alpha = 0.5)
                     
p <- p + geom_hline(yintercept=log(3), lty = 2)                    
p <- p + facet_grid(. ~ fHabitat, scales = "fixed")
p

















#End of exercise
######################################################
#Code below is slightly older code...and used lattice.
#It will not be discussed during the course
#Creat a grid with covariate values
ND <- expand.grid(LCC = seq(0.67, 48.35, length = 10),
                  fHabitat = levels(CO$fHabitat))

#Predict does not work..so manually calculate 
#predicted values and SEs
X  <- model.matrix(~ LCC + fHabitat, data = ND)
mu <- X %*% fixef(M2)
SE <- sqrt(diag(X %*% vcov(M2) %*% t(X)))

#Glue them all together
ND <- cbind(ND, mu, SE)
head(ND)

par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = CO$LCC,
     y = CO$LogG,
     xlab = "LCC",
     ylab = MyYlab,
     pch = 16,
     cex.lab = 1.5)
abline(h = log(3), lty = 2)

for (i in levels(CO$fHabitat)) { 
  NDi <- ND[ND$fHabitat == i, ]	
  lines(NDi$LCC, NDi$mu, lwd = 3) 
  }   



#Or very advanced plotting code....but very (!) nice for a paper
#Get predicted values and store them in a vector

#Run each time from here
Pred    <- NULL
SE      <- NULL
ID      <- NULL
LCC     <- NULL
for (i in levels(CO$fHabitat)){
   MinMaxi <- range(CO$LCC[CO$fHabitat==i])
   MyDatai <- data.frame(LCC=
                  seq(from = MinMaxi[1],
                      to   = MinMaxi[2],
                      length=100),
                  fHabitat = factor(i,levels(CO$fHabitat)) )
   Xi <- model.matrix(~ LCC + fHabitat, data = MyDatai)               
   Pi <- Xi %*% fixef(M2)
   Se <- sqrt(diag(Xi %*% vcov(M2) %*% t(Xi)))
   Pred <- c(Pred, Pi)
   SE   <- c(SE, Se)
   ID   <- c(ID, rep(i, nrow(MyDatai)))
   LCC <- c(LCC, MyDatai$LCC)
}

data.frame(LCC, Pred, SE, ID)


xyplot(LogG ~ LCC | fHabitat,
       data = CO,
       xlab = list("LCC", cex = 1.5),
       ylab = list(MyYlab, cex = 1.5),
       layout = c(5,1),   #Modify
       strip = function(bg = 'white', ...)
       strip.default(bg = 'white', ...),
       scales = list(alternating = FALSE,
                     x = list(relation = "same"),
                     y = list(relation = "same")),
       panel=function(x, y, subscripts){
             WhichPanel <- CO$fHabitat[subscripts][1]
             xi  <- LCC[ID == WhichPanel]
             yi  <- Pred[ID == WhichPanel]
             sei <- SE[ID == WhichPanel]
             panel.grid(h=-1, v= 2)
             #panel.lines(xi, yi + 1.96 * sei, col = 1, lty = 2)
             #panel.lines(xi, yi - 1.96 * sei, col = 1, lty = 2)        
             panel.polygon(c(xi, rev(xi)),
                           c(yi - 1.96 * sei, rev(yi + 1.96 * sei)),
                           col = grey(0.5), border = NULL,
                           density = 50)

             panel.lines(xi, yi, col = 1)
             panel.points(x, y, col = 1, pch = 16)
             panel.abline(h=log(3),lty = 2)
             })
#Nice graph...but not easy to make!       










