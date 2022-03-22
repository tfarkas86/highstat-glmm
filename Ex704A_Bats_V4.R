
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

M1 <- glmer(LACI ~ TempMax.std + fDistRoad + fLight100  +
                       fSite + fYear +
                       fSite : fDistRoad + fDistRoad : TempMax.std + 
                       (1 | fTransect) + (1 | fDate),
            family = poisson,
            data = Bats)

summary(M1) #Seems the random effects do their job.
drop1(M1, test = "Chi")


#Check for overdispersion
E1 <- resid(M1, type = "pearson")
N  <- nrow(Bats)
p  <- length(fixef(M1)) + 2
sum(E1^2) / (N - p)
#Overdispersed!


#Why is there overdispersion? Apply model validation
plot(y = E1,
     x = Bats$TempMax.std)
abline(h = 0)

plot(y = E1,
     x = Bats$WindMax)
abline(h = 0)

xyplot(E1 ~ TempMax.std | fTransect,
       data = Bats,
       panel = function(x,y) {
       	panel.points(x,y, pch = 16, col = 1)
       	mi <- lm(y~x)
       	panel.abline(mi, lwd = 2, col = 1)
       })
#Scope for a random slope?
#We can't do the same for Date as on the temp is always the same for all observations from the same day

boxplot(E1 ~ fDistRoad, data = Bats)
boxplot(E1 ~ fLight100, data = Bats)
boxplot(E1 ~ fSite, data = Bats)
boxplot(E1 ~ fYear, data = Bats)

#So...no clear cause for the overdispersion.
#Solutions:
#1. Observation level random intercept
#2A. NB GLMM via nb.glmer
#2B. NB GLMM via glmmADMB
#2C. NB GLMM via JAGS
#3.  Zero inflated Poisson GLMM


##############################################################
#Solution 1. Add observation level random intercept
Bats$Eps <- 1:nrow(Bats)
M2 <- glmer(LACI ~ TempMax.std + fDistRoad + fLight100  +
                       fSite + fYear +
                       fSite : fDistRoad + fDistRoad : TempMax.std + 
                       (1 | fTransect) + (1 | fDate) + (1|Eps),
            family = poisson,
            data = Bats)

drop1(M2, test = "Chi")
summary(M2)





######################################################################
#Solultion 2A. NB GLMM using glmer.nb or glmmadmb
#This approach gives trouble in the sense that some
#standard errors are identical!

#Solution 2B: NB GLMM using glmmadmb
library(glmmADMB)
M3 <- glmmadmb(LACI ~ 1 + TempMax.std + fDistRoad + fLight100  + fSite + fYear + 
                      fSite : fDistRoad + fDistRoad : TempMax.std + 
                      (1 | fTransect) + (1 | fDate),
               data = Bats,
               family = "nbinom")
summary(M3)

#Check overdispersion
E3 <- resid(M3, type = "pearson")
N  <- nrow(Bats)
p  <- length(fixef(M3)) + 2
sum(E3^2) / (N - p)
#Perfect


#Do a model selection on the interactions
#Drop fSite : fDistRoad
M3A <- glmmadmb(LACI ~ 1 + TempMax.std + fDistRoad + fLight100  + fSite + fYear + 
                      fDistRoad : TempMax.std + 
                      (1 | fTransect) + (1 | fDate),
               data = Bats,
               family = "nbinom")

#Drop fDistRoad : TempMax.std
M3B <- glmmadmb(LACI ~ 1 + TempMax.std + fDistRoad + fLight100  + fSite + fYear + 
                      fSite : fDistRoad  + 
                      (1 | fTransect) + (1 | fDate),
               data = Bats,
               family = "nbinom")

AIC(M3, M3A, M3B)
    # df     AIC
# M3  17 782.300
# M3A 13 786.998
# M3B 15 781.562  <----




#Next round of model selection
M4 <- glmmadmb(LACI ~ 1 + TempMax.std + fDistRoad + fLight100  + fSite + fYear + 
                      fSite : fDistRoad  + 
                      (1 | fTransect) + (1 | fDate),
               data = Bats,
               family = "nbinom")

#Drop fSite : fDistRoad
M4B <- glmmadmb(LACI ~ 1 + TempMax.std + fDistRoad + fLight100  + fSite + fYear  + 
                      (1 | fTransect) + (1 | fDate),
               data = Bats,
               family = "nbinom")

AIC(M4, M4B)
    # df     AIC
# M4  15 781.562 <--
# M4B 11 787.220  

#You could go on with model selection...if you want to.
#We stop here and focus on model M4

#Is everything significant?
drop1(M4, test = "Chi") 
#                Df    AIC    LRT  Pr(>Chi)    
# <none>             781.56                     
# TempMax.std      1 783.09  3.524 0.0604863 .  
# fLight100        1 780.87  1.312 0.2520329    
# fYear            1 793.13 13.570 0.0002298 ***
# fDistRoad:fSite  4 787.22 13.658 0.0084705 ** 

# Well...we now answered a major question...there
# is no light effect. 



summary(M4)




#####################################################################################
#Model validation
E4 <- resid(M4, type = "pearson") 

#Task: plot residuals vs each covariate
#      Plot residuals vs fitted values




#####################################################################################
#Model interpretation
#Task: Write down the fitted model


#This was the model
#M4 <- glmmadmb(LACI ~ 1 + TempMax.std + fDistRoad + fLight100  + fSite + fYear + 
#                      fSite : fDistRoad  + 
#                      (1 | fTransect) + (1 | fDate),
#               data = Bats,
#               family = "nbinom")

# Let's visualize the effect of the two most important components, 
# year and the interaction fSite : fDistRoad.

#Create artificial covariate values
#Average temperature, and no light.
MyData <- expand.grid(TempMax.std = 0,
                      fDistRoad = levels(Bats$fDistRoad),
                      fLight100 = factor("No", levels = c("No", "Yes")),
                      fSite     = levels(Bats$fSite),
                      fYear     = levels(Bats$fYear))

#Convert the covariate values into an X matrix
X4 <- model.matrix(~ TempMax.std + fDistRoad + fLight100  + fSite + fYear + 
                     fSite : fDistRoad,
                     data = MyData) 

#Extract estimated parameters
beta4    <- fixef(M4)
CovBeta4 <- vcov(M4)

#And calculate the fitted values and the 95% CIs
MyData$eta  <- X4 %*% beta4
MyData$Fit  <- exp(X4 %*% beta4)
MyData$SE   <- sqrt( diag (X4 %*% vcov(M4) %*% t(X4) ) ) 
MyData$seup <- exp(MyData$eta + 1.96 * MyData$SE)
MyData$selo <- exp(MyData$eta - 1.96 * MyData$SE)

library(ggplot2)  
p <- ggplot()
p <- p + xlab("Feeding type") + ylab("Shell length")
p <- p + theme(text = element_text(size=15)) + theme_bw()

p <- p + geom_point(data = MyData, 
                    aes(x = fDistRoad, 
                        y = MyData$Fit, 
                        size = 6),    
                    col = ("black"))

p <- p + geom_errorbar(data = MyData,
                       aes(x = fDistRoad, 
                           ymax = seup, 
                           ymin = selo), 
                       width=0.2)

p <- p + geom_point(data = Bats, 
                    aes(x = fDistRoad, y = LACI),
                    position = position_jitter(width = .02),
                    color = grey(0.3),
                    size = 2)
       

                    
p <- p + facet_grid(fYear ~ fSite, 
                    scales = "fixed")
p <- p + theme(legend.position="none") 
p <- p + theme(strip.text.y = element_text(size = 15, 
                                           colour = "black", 
                                           angle = 20),
               strip.text.x = element_text(size = 15, 
                                           colour = "black", 
                                           angle = 0)                            
                                           )
p

#Rey to explain what this means
###########################################################


################################################################################
#What have we learned?
# 1. Keep your models simple. We have a small data set with lots of parameters!



