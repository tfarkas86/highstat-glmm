
#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
# Fungi exercise
# These data are taken from:
# Elucidating the nutritional dynamics of fungi using 
# stable isotopes
# Jordan R. Mayor, Edward A. G.Schuur and Terry W. Henkel
# Ecology Letters, (2009) 12: 171-183


#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/NewDataSets/lme_Fungi/")
Fungi <- read.csv("FungiV5.csv", 
                   header = TRUE,
                   dec = ".")
names(Fungi)
 # [1] "Study"    "Species"  "Site"     "Temp"     "Rain"     "Type"    
 # [7] "d13C"     "d15N"     "Category" "STND_13c" "STND_15n" "Lat"     
# [13] "Long"    

str(Fungi)

# Fungi function at two fundamental biogeochemical 
# interfaces between soil and plants. 
#   -Decomposer fungi mineralize organic carbon compounds 
#    in detritus and liberate mineral nutrients in the process.
#   -Mycorrhiza-forming fungi, as mutualistic root extensions, 
#    enhance mineral, and perhaps organic, nutrient uptake in 
#    exchange for plant photosynthate

# Fungi are divided into saprotrophic (SAP) and ectomycorrhizal (ECM) 
# functional groups

# Mycorrhizal and saprotrophic (SAP) fungi are essential 
# to terrestrial element cycling due to their uptake of 
# mineral nutrients and decomposition of detritus. 

# d15N and d13C of fungi provide time integrated biogeochemical# information regarding the acquisition, transformation,# and export of C and N by fungi.

# In order to determine if fungal d15N or d13C# patterns across ecosystems are similar to plants and 
# soils, we assessed the explanatory capacity of:
#   -mean annual temperature (Temp...called MAT in paper),
#   -mean annual precipitation (Rain..called MAP in paper), and 
#   -latitude (Lat). 

#Methods:
# 1. The athors compiled d15N and d13C values from 
#    one novel and ten published data sets. Compiled 
#    data included 913 d15N and 813 d13C values from 
#    collector categorized ECM or SAP fungi, and 27 
#    fungi of unknown ecological role, together comprising 
#    148 genera.
#2.  32 study sites
###################################



###########################################################
#Load packages and library files
library(lattice)  #Needed for multi-panel graphs
library(nlme)
library(mgcv)
library(ggplot2)

#source the file: HighstatLibV7.R
source("/Users/Highstat/MaggieWindows/applicat/HighlandStatistics/webHighstatNew2_2010/CourseMCMCGLMM/RCode/HighstatLibV7.R")
##################################################################


# Task:
# Model d13C and d15N as a function of:
# Temp, Rain, Lat, Type (ECM, SAP, Unknown) + relevant interactions
# Random effect site 
###################################################




###################################################
#Housekeeping
Fungi$fSite    <- factor(Fungi$Site)

###################################################



###################################################
#Data exploration
#Outliers
MyVar <- c("Temp", "Rain", "Lat", "Long", "d15N", "d13C")
Mydotplot(Fungi[,MyVar])
# One small d13C value
# Lots of observations with the same covariate values!

#Are the data balanced?
table(Fungi$fSite)
sort(table(Fungi$Species))
table(Fungi$Study)
table(Fungi$Type)


#Collinearity
boxplot(Temp ~ Type, 
        data = Fungi, 
        xlab = "Type", 
        ylab = "Temp")
#Small collinearity?

boxplot(Rain ~ Type, 
        data = Fungi, 
        xlab = "Type", 
        ylab = "Rain")
#Minor collinearity?

#Site will be collinear with temp and rain.
#And latitude


MyVar <- c("d13C", "d15N", "Temp", "Rain","Lat")
Mypairs(Fungi[, MyVar])
#Serious collinearity! Pick either Temp, rain or lat!
#Pick Rain?

p <- ggplot()
p <- p + geom_point(data = Fungi, 
                    aes(y = d15N, x = Rain),
                    shape = 16, 
                    size = 3)
p <- p + xlab("Rain") + ylab("d15N")
p <- p + theme(text = element_text(size=15))
p <- p + geom_smooth(data = Fungi, method = "lm", 
                   aes(x = Rain, y = d15N), 
                   colour = "black")                     
p <- p + facet_grid(. ~ Type, scales = "fixed")
p

#There is probably an interaction. But the quality of the
#rain data in the unknown group is not that nice!
#Should we analyse the data without the unknowns?


####################
#Relationships
boxplot(d15N ~ fSite, 
        data = Fungi, 
        xlab = "Site")
#Site effect!
        
boxplot(d15N ~ Study, 
        data = Fungi, 
        xlab = "Study") 
#Location effects?

boxplot(d15N ~ Species, 
        data = Fungi, 
        xlab = "Species")
#Not many observations per species


par(mar = c(5,6,2,2), cex.lab = 1.5)
plot(x = Fungi$Rain, 
     y = Fungi$d15N,
     pch = 16,
     cex = 1,
     xlab = "Rain",
     ylab = "d15N")

par(mar = c(5,6,2,2), cex.lab = 1.5)
plot(x = Fungi$Temp, 
     y = Fungi$d15N,
     pch = 16,
     cex = 1,
     xlab = "Temperature",
     ylab = "d15N")


par(mar = c(5,6,2,2), cex.lab = 1.5)
plot(x = Fungi$Lat, 
     y = Fungi$d15N,
     pch = 16,
     cex = 1,
     xlab = "Latitude",
     ylab = "d15N")
     
     
     
#Dump NAs
MyVar <- c("d13C", "d15N", "Temp", "Rain", 
           "Lat", "Long", "fSite", "Type",
           "Study", "Species")
Fungi2 <- na.exclude(Fungi[, MyVar])
dim(Fungi)
dim(Fungi2)


table(Fungi2$fSite, Fungi2$Type)
#it would make sense to do only the ECM and SAP data!
#############


#Data exploration indicates that:
# There are is potentially 1 outlier
# There is serious collinearity
# There could be an interaction
# Should we drop the Unknowns?


##########################################################



#Start of analysis. First the frequentist analysis

M0 <- lme(d13C ~ Rain *  Type,
          random =~ 1 | fSite,
          data = Fungi2, 
          method = "REML")
summary(M0)
#By design of the study we use site as random effect.
             
                   
                       
#Model validation                       
E0 <- resid(M0, type = "n")
F0 <- fitted(M0)

par(mfrow = c(3,3), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = F0, 
     y = E0, 
     xlab = "Fitted values",
     ylab = "Normalized residuals")
abline(h = 0, lty = 2)

plot(x = Fungi2$Temp, 
     y = E0, 
     xlab = "Temperature",
     ylab = "Normalized residuals")
abline(h = 0, lty = 2)

plot(x = Fungi2$Rain, 
     y = E0, 
     xlab = "Rainfall",
     ylab = "Normalized residuals")
abline(h = 0, lty = 2)

plot(x = Fungi2$Lat, 
     y = E0, 
     xlab = "Lat",
     ylab = "Normalized residuals")
abline(h = 0, lty = 2)


boxplot(E0 ~ Species,   
        ylab = "Normalized residuals",
        data = Fungi2, 
        xlab = "Species")
abline(h = 0, lty = 2)

boxplot(E0 ~ Study, 
        ylab = "Normalized residuals", 
        data = Fungi2,  
        xlab = "Study")
abline(h = 0, lty = 2)

boxplot(E0 ~ fSite, 
         ylab = "Normalized residuals",
         data = Fungi2, 
         xlab = "Site",
         varwidth = TRUE)
abline(h = 0, lty = 2)


plot(x = F0,
     y = Fungi2$d13C,
     xlab = "Fitted values",
     ylab = "Residuals")
#####################################################          


#Step 2: Are all terms significant?
M1 <- lme(d13C ~ Rain *  Type,
          random =~ 1 | fSite,
          data = Fungi2, 
          method = "ML")
summary(M1)
M1A <- update(M1, .~. - Rain :  Type)
anova(M1, M1A)
#Yes....significant.



#Step 3: Explain what it all means.
summary(M0)  #This is REML obtained stuff
#Calculate the intra-class correlation.
0.8901247 ^2 / (0.8901247 ^2 +  1.214055 ^2)



#TASK: Write down the 3 equations




##############################################
#Sketch fitted values
#A. Specify covariate values for predictions
#B. Create X matrix with expand.grid
#C. Calculate predicted values
#D. Calculate standard errors (SE) for predicted values
#E. Plot predicted values
#F. Plot predicted values +/- 	1.96 * SE

library(plyr)
MyData <- ddply(Fungi2, .(Type), summarize,
                Rain = seq(min(Rain), max(Rain),length = 25))
MyData


X  <- model.matrix(~ Rain *  Type, data = MyData)
mu <- X %*% fixef(M0)
SE <- sqrt(diag(X %*% vcov(M0) %*% t(X)))


#Glue them all together
MyData$mu <- mu
MyData$ub <- mu + 1.96 * SE
MyData$lb <- mu - 1.96 * SE
head(MyData)


p <- ggplot()
p <- p + geom_point(data = Fungi, 
                    aes(y = d13C, x = Rain),
                    shape = 16, 
                    size = 3)
p <- p + xlab("Rain") + ylab("d13C")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = Rain, y = mu), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = Rain, 
                         ymax = ub, 
                         ymin = lb ),
                     alpha = 0.5)
p <- p + facet_grid(. ~ Type, scales = "fixed")
p


#What happens if we drop the unknowns?
Fungi3 <- Fungi2[Fungi2$Type != "Unknown", ]
dim(Fungi3)
FungiBackup <- Fungi2
Fungi2 <- Fungi3  #Now rerun the code
Fungi2$Type <- factor(Fungi2$Type)


