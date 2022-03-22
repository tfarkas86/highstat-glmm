#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
#Cowbird exercise
#################################################################
#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/FilesMCMC_GLMM/CDMCMCGLMGLMMCourse/Data")
CB <- read.table(file = "CowbirdV3Book.txt",
                 header = TRUE)
dim(CB)
names(CB)
str(CB)
#################################################################

#Task: 
#Some bird species place their eggs solely in the nests 
#of other bird species; this is called brood parasitism. 
#Brood parasitism has led to an evolutionary battle between 
#host and parasite; the parasite tries to fool the host 
#(e.g., mimic the host eggs), and the host tries to outsmart 
#the parasite with advanced defensive mechanisms (e.g., nest location). 
#Patten et al. (2011) used an extensive dataset and looked for 
#cues in brood parasitism nest selection by the brood parasitic 
#brown-headed Cowbird (Molothrus ater).

#The dataset consists of a large number of observations (2,825) 
#on absence and presence of brood parasitism by cowbirds. 
#Patten et al. (2011) used logistic regression to investigate 
#whether the probability of brood parasitism is affected by 
#covariates such as perch proximity, nest height, livestock 
#proximity, habitat density, nest exposure, laying date of 
#the first egg, host clutch size, and host species, among other 
#variables. They found egg laying date to have a significant effect.
#The question we focus on in this chapter builds on the work 
#by Patten and colleagues: 

#Is the effect of egg laying date on the probability of brood 
#parasitism by cowbirds the same for all host species, or does this 
#effect differ per host species? 
#And if the effect differs per species, how does it differ from 
#species to species?

###################################################################
#Load packages and library files
library(lattice)  #Needed for multi-panel graphs
library(lme4)
source(file = "/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV9.R")  
########################################################




########################################################
#Data exploration

#First we dump a few NAs
CB1 <- na.exclude(CB)
dim(CB1)

#Outliers
dotchart(CB1$Firstegg)
par(mar = c(5,5,4,2))
plot(x = CB1$Firstegg,
     y = 1:nrow(CB1),
     pch = 16,
     cex = 0.5,
     xlab = "Range of the data",
     ylab = "Order of the data",
     cex.lab = 1.5,
     main = "Less memory in Word")
#That is ok

dotchart(CB1$Cowbird)


table(CB1$Cowbird)
table(CB1$Plot)
table(CB1$Year)
#That is all ok-ish

table(CB1$Species)
#This is not ok!! If we want to include interactions, then
#we should dump all species < 25 or 30-ish obervations

#These data were used for a GAMM chapter, and we included
#a Firstegg x Species smoother interaction. For this 
#you need at least 50-ish or more number
#of observations per species.

#This is what we used in the 2014 book
CB2 <- subset(CB1, Species == "DICK" |
                   Species == "EAME" |
                   Species == "GRSP" |
                   Species == "RWBL")
CB2$Species <- factor(CB2$Species)
#But for GLMM techniques you could as well
#add COGR, EAPH and OROR
#Using only 4 species keeps computing time
#and coding of GAMM relatively low and simple.


#Check whether this is still ok:
table(CB2$Cowbird)
table(CB2$Plot)
table(CB2$Year)
#Yes...that is all ok

     
#Collinearity
boxplot(Firstegg ~ Year, data = CB2) 
boxplot(Firstegg ~ Species, data = CB2) 
boxplot(Firstegg ~ Plot, data = CB2) 
#Some minor collinearity going on!

    
#Relationships    
plot(x = CB2$Firstegg,
     y = CB2$Cowbird)
#That is not very useful.
#Try this:
boxplot(Firstegg ~ Cowbird, data = CB2)     
boxplot(Firstegg ~ Cowbird * Species, data = CB2)     



with(CB2, table(Species, Cowbird))
with(CB2, table(Year, Cowbird))
with(CB2, table(Species, Year, Cowbird))


           
#Balanced data?     
tapply(CB2$Firstegg, INDEX = CB2$Species, FUN = range)
# $DICK
# [1]  45 133
# $EAME
# [1]  19 114
# $GRSP
# [1]  34 131
# $RWBL
# [1]  26 116

#That could be trouble!
#Since we are interested in the interaction
#between species and Firstegg, it may be an option
#to focus the analysis on on the Firstegg data between
#days 45 and 114
#CB3 <- subset(CB2, Firstegg >= 45 & Firstegg <= 114)

#Alternatively..use the data as it is...and be careful 
#with the interpretation. That's what we will do here
CB3 <- CB2


################################################
#Some housekeeping
CB3$fSpecies <- factor(CB3$Species)
CB3$fYear    <- factor(CB3$Year)
CB3$fPlot    <- factor(CB3$Plot)


#######################################################
#By design of the study use fPlot as random effect
#Step 1: Apply a binomial GLMM
#detach(package:nlme)
CB3$Firstegg.std <- (CB3$Firstegg - mean(CB3$Firstegg)) / sd(CB3$Firstegg)

M1 <- glmer(Cowbird ~ Firstegg * fSpecies + fYear + (1|fPlot),
            data = CB3, 
            family = binomial)

M1b <- glmer(Cowbird ~ Firstegg * fSpecies + fYear + (1|fPlot),
            data = CB3, 
            family = binomial,
            control = glmerControl(optimizer="bobyqa"))

M1c <- glmer(Cowbird ~ Firstegg * fSpecies + fYear + (1|fPlot),
            data = CB3, 
            family = binomial,
            control = glmerControl(optimizer="Nelder_Mead"))


summary(M1)
summary(M1b)
summary(M1c)

drop1(M1, test = "Chi")

# Fixed effects:
                       # Estimate Std. Error z value Pr(>|z|)   
# (Intercept)           -0.593780   0.435904  -1.362  0.17314   
# Firstegg              -0.012481   0.004448  -2.806  0.00502 **
# fSpeciesEAME          -0.221750   0.642690  -0.345  0.73007   
# fSpeciesGRSP          -0.605255   0.861386  -0.703  0.48227   
# fSpeciesRWBL           0.048746   0.597494   0.082  0.93498   
# fYear1993              0.689879   0.234533   2.942  0.00327 **
# fYear1994              0.157182   0.235826   0.666  0.50508   
# fYear1995              0.515807   0.233966   2.205  0.02748 * 
# fYear1996              0.437915   0.229227   1.910  0.05608 . 
# Firstegg:fSpeciesEAME -0.030859   0.012293  -2.510  0.01206 * 
# Firstegg:fSpeciesGRSP -0.011287   0.012822  -0.880  0.37870   
# Firstegg:fSpeciesRWBL -0.009959   0.008375  -1.189  0.23441   


# Cowbird_ij ~ B(p_ij)
# E(Cowbird_ij) = p_ij

#           exp( eta_ij )
# p_ij = -------------------
#          1 + exp( eta_ij )

# eta_ij = Fixed stuff + a_i

# What is eta_ij for DICK & 1992? 
# eta_ij = -0.593 - 0.012 * FirstEgg_ij + a_i

# What is eta_ij for RWBL & 1996? 
# eta_ij = -0.593 + 0.048 + 0.437915 +  
#         (-0.012 - 0.0099) * FirstEgg_ij + a_i
###################################################


#R2 for all basic glmm models (gaussian, poisson, binomial)
# marginal is for fixed effects only, conditional for fixed + random
source(file = "R2glmms.R")
r.squared.merMod(M1) 
     # Class   Family  Link Marginal Conditional      AIC
# 1 glmerMod binomial logit 0.158565   0.3557762 1916.253



###################################################
#Model validation of binomial GLMs and GLMMs is difficult 	
E1 <- resid(M1, type = "pearson")
plot(x = CB3$Firstegg,
     y = E1)

#Step 2: Everything is significant

#Step 3: Explain what it all means.

#Create data on grid, and the matching X matrix
range(CB3$Firstegg)
MyData <- expand.grid(fYear    = levels(CB3$fYear),
                      Firstegg = seq(19, 133, length = 25),
                      fSpecies = levels(CB3$fSpecies))
X <- model.matrix(~Firstegg * fSpecies + fYear, data = MyData)
#Extract parameters and parameter covariance matrix
betas    <- fixef(M1)
Covbetas <- vcov(M1)

#Calculate the fitted values in the predictor scale
MyData$eta <- X %*% betas
MyData$Pi  <- exp(MyData$eta) / (1 + exp(MyData$eta))

#Calculate the SEs on the scale of the predictor function
MyData$se    <- sqrt(diag(X %*% Covbetas %*% t(X)))
MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$se) / (1 + exp(MyData$eta  + 1.96 *MyData$se))
MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$se) / (1 + exp(MyData$eta  - 1.96 *MyData$se))

head(MyData)



library(ggplot2)
p <- ggplot()
p <- p + geom_point(data = CB3, 
                    aes(y = Cowbird, x = Firstegg),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Date of first egg") + ylab("Probability of parasitism")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = Firstegg, y = Pi), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = Firstegg, 
                         ymax = SeUp, 
                         ymin = SeLo ),
                     alpha = 0.2)
p <- p + facet_grid(fYear ~ fSpecies, scales = "fixed")
p


#And now multiple species in the same panel
p <- ggplot()
p <- p + geom_point(data = CB3, 
                    aes(y = Cowbird, x = Firstegg),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Date of first egg") + ylab("Probability of parasitism")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = Firstegg, 
                       y = Pi, 
                       group = fSpecies,
                       colour = fSpecies))

p <- p + geom_ribbon(data = MyData, 
                     aes(x = Firstegg, 
                         ymax = SeUp, 
                         ymin = SeLo,
                         group = fSpecies,
                         fill = fSpecies),
                         alpha = 0.2)
p <- p + facet_grid(. ~ fYear, scales = "fixed")
p
###########################################################


#And now multiple years in the same panel
p <- ggplot()
p <- p + geom_point(data = CB3, 
                    aes(y = Cowbird, x = Firstegg),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Date of first egg") + ylab("Probability of parasitism")
p <- p + theme(text = element_text(size=15)) + theme_bw()
p <- p + geom_line(data = MyData, 
                   aes(x = Firstegg, 
                       y = Pi, 
                       group = fYear,
                       colour = fYear))

p <- p + geom_ribbon(data = MyData, 
                     aes(x = Firstegg, 
                         ymax = SeUp, 
                         ymin = SeLo,
                         group = fYear,
                         fill = fYear),
                         alpha = 0.2)
p <- p + facet_grid(. ~ fSpecies, scales = "fixed")
p
###########################################################







###########################################################
#And the lattice equivalent
xyplot(Cowbird ~ Firstegg | fYear,
       data = CB3,
       xlab = "First egg",
       ylab = "Probability of parasitism",
       panel = function(x,y,subscripts,...) {
         panel.points(x, y, col = 1, pch = 16)
       	 ThisPanel <- CB3$fYear[subscripts][1]
         MyDatai <- MyData[MyData$fYear == ThisPanel, ]
         AllSpecies <- levels(MyData$fSpecies)
         for (i in AllSpecies) {
        		xi <- MyDatai[MyDatai$fSpecies == i, "Firstegg"]
        		yi <- MyDatai[MyDatai$fSpecies == i, "Pi"]
        		panel.lines(xi, yi, col = 1)
       	}})



