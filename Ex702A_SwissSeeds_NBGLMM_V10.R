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
length(table(SPS$Origin))
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
# BMtL10mg	vegetative growth		mg
# BMS10mg  	reproductive biomass		mg
# H10f	    number of inflorescences	
# S10f	    number of flower stalks	
# tBB10d	flowering initiation	  	DOI
# tSB10d	first mature fruit	  	DOI
# dB10d	    flowering duration		days

#Due to vole damage to the plants, one site 
#needs to be dropped from the analysis. This is:
#Neus

SPS2 <- SPS[SPS$Site != "Neus",]
dim(SPS)
dim(SPS2)



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
library(lattice)
source("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV7.R")
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
SPS3        <- SPS2[!is.na(SPS2$H10f) & !is.na(SPS2$InitSize) , MyVar]
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
tapply(SPS3$InitSize, FUN = range, INDEX = SPS3$Alt)
tapply(SPS3$InitSize, FUN = range, INDEX = SPS3$Soil)

#This is a sensible option:
SPS4 <- SPS3[SPS3$InitSize <= 2.2,]
dim(SPS3)
dim(SPS4)
#################################################





#################################################
#Start frequentist analysis
#Fit the following model
# H10f is the number of inflorescences
# H10f_ijk ~ Poisson(mu_ijk)
# E(H10f_ijk) = mu_ijk
# log(mu_ijk) = InitSize_ijk + 
#               Soil treatment_ijk * altitude treatment_ijk + 
#               Origin effect_l + Site effect_i + 
#               Bed-within-site effect_ij


library(lme4)


#Apply Poisson GLMM
SPS4$InitSizec <- (SPS4$InitSize - mean(SPS4$InitSize)) / sd(SPS4$InitSize)
M1 <- glmer(H10f ~ InitSizec + Soil * Alt +
                  (1| Origin) + (1 | Site / Bed), 
                  data = SPS4,
                  family = poisson)
summary(M1)
E1 <- resid(M1, type = "pearson")
N  <- nrow(SPS3)
p  <- length(fixef(M1)) + 3
sum(E1^2) / (N - p)
#Small overdispersion


#Apply NB GLMM
M2 <- glmer.nb(H10f ~ InitSizec + Soil * Alt +
                  (1| Origin) + (1 | Site/Bed), 
                  control = glmerControl(optimizer="bobyqa"),
                  data = SPS4)
summary(M2)
E2 <- resid(M2, type = "pearson")
N  <- nrow(SPS3)
p  <- length(fixef(M2)) + 3 + 1
sum(E2^2) / (N - p)

#Warning???


#Check whether glmmADMB gives the same results
library(glmmADMB)
M3 <- glmmadmb(H10f ~ InitSizec + Soil * Alt +
                      (1| Origin) + (1 | Site/Bed),
               family = "nbinom",       
               data = SPS4)

E3 <- resid(M3, type = "pearson")
N  <- nrow(SPS4)
p  <- length(fixef(M3)) + 3 + 1
sum(E3^2) / (N - p)

summary(M3)
summary(M2)

#We can also plot them with coefplot2
library(coefplot2)
coefplot2(M2, intercept = TRUE)
coefplot2(M3, intercept = TRUE, add = TRUE, offset = 0.1)

#So...both packages give similar results.



# Task: Try to write down the model that we just fitted.






AIC(M1, M2, M3)
   # df      AIC
# M1 10 3998.640 Poisson GLMM
# M2 11 3121.632 NB GLMM with glmer.nb
# M3 11 3121.700 NB GLMM with glmmadmb


#Now drop the interaction from M3
M3a <- glmmadmb(H10f ~ InitSizec + Soil + Alt +
                      (1| Origin) + (1 | Site/Bed),
               family = "nbinom",       
               data = SPS4)


AIC(M3, M3a) #We might as well drop the interaction!
summary(M3a)

#Should we do a model selection on the main terms as well,
#or keep it as it is?
#If you want to know whether these main terms are significant,
#then drop them, get the likelihood and manually calculate a
#p-value. Or hope that drop1 works correctly for a ZINB GLMM


#Homework: Apply model validation again


####################################################################
#Sketch fitted values
#A. Specify covariate values for predictions
#B. Create X matrix with expand.grid
#C. Calculate predicted values
#D. Calculate standard errors (SE) for predicted values
#E. Plot predicted values
#F. Plot predicted values +/- 2* SE 


#A: 
range(SPS4$InitSizec)
#0.2  2.2
MyData <- expand.grid(InitSizec = seq(-1.84, 2.77, length = 25),
                      Soil      = levels(SPS4$Soil),
                      Alt       = levels(SPS4$Alt))
                       
#B. Create X matrix with expand.grid
X <- model.matrix(~ InitSizec + Soil + Alt, data = MyData)

#C. Calculate predicted values
MyData$eta  <- X %*% fixef(M3a) 
MyData$ExpY <- exp(X %*% fixef(M3a))


#D. Calculate standard errors (SE) for predicted values
#SE of fitted values are given by the square root of
#the diagonal elements of: X * cov(betas) * t(X)
MyData$SE <- sqrt(  diag(X %*%vcov(M3a) %*% t(X))  )

MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$SE) 
MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$SE) 




library(ggplot2)
p <- ggplot()
p <- p + geom_point(data = SPS4, 
                    aes(y = H10f, x = InitSizec),
                    shape = 1, 
                    size = 1)
p <- p + xlab("InitSizec") + ylab("H10f")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = InitSizec, y = ExpY), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = InitSizec, 
                         ymax = SeUp, 
                         ymin = SeLo ),
                     alpha = 0.2)
p <- p + facet_grid(Soil ~ Alt, scales = "fixed")
p

#Note the points with the poor fit in the middle.
#I smoother of InitSizec would be able to capture it.
#I would try a ZIP GAM...where the smoother can capture
#the non-linear pattern of InitSize and also a InitiSize
#in the binary part to capture the zeros on the right hand side
#of the gradient? #But that requires some fancy coding.

