#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
#Coral reef exercise
#These data are taken from:
#Caribbean-wide decline in carbonate production threatens 
#coral reef growth#Perry et al. (2012)
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
source("/Users/Highstat/MaggieWindows/applicat/HighlandStatistics/webHighstatNew2_2010/CourseMCMCGLMM/RCode/HighstatLibV7.R")
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


# The response variable used in Perry et al. is 
# Net_production, which is calculated as:

# Net_production  = Gross_production - Gross_erosion

#Instead of Net_production we can also look at:

# Ratio = Gross_production / Gross_erosion
 
# A value > 1 means an increase and a value <1 means a decrease
#
# We have one observation within a transect. 
# Hence, transect is the statistical sampling unit. 
# Potential explanatory variables are Country, 
# Site (a group of transects), Depth, Habitat, and 
# LCC (life coral cover %). 

# Task:
# Model the ratio as a function of:
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
CO$Ratio <- CO$Gross_production / CO$Gross_erosion
#CO$G     <- CO$Net_production

###################################################



###################################################
#Data exploration
#Outliers
MyVar <- c("Ratio", "LCC")
Mydotplot(CO[,MyVar])


#Are the data balanced?
table(CO$fCountry)
table(CO$fSite)
table(CO$fHabitat)

#4 countries
#14 uniquely defined sites


#Relationships
boxplot(Ratio ~ fSite, data = CO, xlab = "Site")
boxplot(Ratio ~ fCountry, data = CO, xlab = "Country")
boxplot(Ratio ~ fHabitat, data = CO, xlab = "Habitat")

MyYlab <- "Ratio of gross production and erision"
par(mar = c(5,6,2,2), cex.lab = 1.5)
plot(x = CO$LCC, 
     y = CO$Ratio,
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
library(lme4)
M0 <- glmer(Ratio ~ LCC + fHabitat + (1 |fSite),
          data = CO,
          family = Gamma(link = "log"))
summary(M0)

#By design of the study we use site as random effect.
#There are only 4 countries, so we can't use a 2-way nested structure.
                   
                       
#Model validation                       
par(mfrow = c(2,2), mar = c(5,5,2,2), cex.lab = 1.4)
E0 <- resid(M0, type = "pearson")
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

#Ouch....do we see a non-linear effect of LCC?
#And heterogeneity?

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


#############################################################




#Step 2: Are all terms significant?
drop1(M0, test = "Chi")
#Everything significant (though we don't know yet whether
#the habitat effect is in fact a country effect)

#Step 3: Explain what it all means.
summary(M0)  

#Write down the 5 equations for Ratio (one per habitat type).


##############################################
#Sketch fitted values
#A. Specify covariate values for predictions
#B. Create X matrix with expand.grid
#C. Calculate predicted values
#D. Calculate standard errors (SE) for predicted values
#E. Plot predicted values
#F. Plot predicted values +/- 	1.96 * SE


#Here is some fancy code to make a better MyData
library(plyr)
MyData <- ddply(CO, .(fHabitat), summarize,
                LCC = seq(min(LCC), max(LCC),length = 10))
MyData



X  <- model.matrix(~ LCC + fHabitat, data = MyData)
eta <- X %*% fixef(M0)
SE <- sqrt(diag(X %*% vcov(M0) %*% t(X)))


#Glue them all together
MyData$mu <- exp(eta)
MyData$ub <- exp(eta + 1.96 * SE)
MyData$lb <- exp(eta - 1.96 * SE)
head(MyData)


library(ggplot2)
p <- ggplot()
p <- p + geom_point(data = CO, 
                    aes(y = Ratio, x = LCC),
                    shape = 16, 
                    size = 3)
p <- p + xlab("LCC") + ylab("Ratio")
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
