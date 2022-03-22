#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#######################################################################

#These data were taken from:
# A resource-based conservation approach for an 
# endangered ecotone species: the Ilex 
# Hairstreak (Satyrium ilicis) in 
# Flanders (north Belgium)

# Dirk Maes, Ilf Jacobs, Natascha Segers, Wouter Vanreusel 
# Toon Van Daele, Guy Laurijssens, Hans Van Dyck

########################################################
#Load the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data")
BF <- read.table("ButterfliesNEggs_V3.txt", header = TRUE)
names(BF)
str(BF)
###########################################################################



###########################################################################
#Load support files and packages
source(file = "HighstatLibV9.R")
library(lme4)
library(lattice)
###########################################################################


# Response: Number of eggs.
# Covariates: We have a whole bunch of covariates, 
#             and their names are self-explanatory.
# We have multiple observation per 'landscape'. That 
# is like a site.

# The variables "Distance2Edge" and "NearestTallOak" are
# expressed in meters.
###########################################################################




###########################################################################
#Data exploration

#Outliers
MyVar <- c("NEggs", "LandscapeID", "NLowBranches",
          "TreeHeight", "HerbCover", "Shelter", 
          "BuckthornAbundance", "BrambleAbundance", "SmallOakAbundance",   
          "Distance2Edge", "NearestTallOak", "X",  "Y" )
Mydotplot(BF[,MyVar])
#Nothing too alien


#Zero inflation!
table(BF[,"NEggs"])


###################################
#Define factors as factors
BF$fShelter     <- factor(BF$Shelter)
BF$fLandscapeID <- factor(BF$LandscapeID)

#How many observations per cluster, and per shelter level?
table(BF$LandscapeID)
table(BF$fShelter)

#Is there a landscape effect?
boxplot(NEggs ~ LandscapeID, data = BF)


#Collinearity
MyVar <- c("NLowBranches",
          "TreeHeight", "HerbCover",
          "BuckthornAbundance", "BrambleAbundance", "SmallOakAbundance",   
          "Distance2Edge", "NearestTallOak", "X",  "Y" )

corvif(BF[,MyVar])
Mypairs(BF[,MyVar])


#Spatial position of the sampling locations
#X and Y are in Belgium coordinates
xyplot(Y ~ X,
       data = BF,
       aspect = "iso")
###################################################



###################################################
#Frequentist analysis

#Standardize the continuous covariates
MyStd <- function(x) { (x - mean(x)) / sd(x)}

BF$NLowBranchesc  <- MyStd(BF$NLowBranches)
BF$TreeHeightc    <- MyStd(BF$TreeHeight)
BF$HerbCoverc     <- MyStd(BF$HerbCover)
BF$BuckthornAbundancec <- MyStd(BF$BuckthornAbundance)
BF$BrambleAbundancec   <- MyStd(BF$BrambleAbundance)
BF$SmallOakAbundancec  <- MyStd(BF$SmallOakAbundance)
BF$Distance2Edgec      <- MyStd(BF$Distance2Edge)
BF$NearestTallOakc     <- MyStd(BF$NearestTallOak)
                    
                    
#Apply a Poisson GLMM using lme4
M1 <- glmer(NEggs ~ NLowBranchesc + TreeHeightc + HerbCoverc + 
                    BuckthornAbundancec + BrambleAbundancec + 
                    SmallOakAbundancec + Distance2Edgec + NearestTallOakc +
                    fShelter +
                    (1 | fLandscapeID),
          family = poisson,
          data = BF)
summary(M1)
#Worrying error messages!

#Check for overdispersion
E1 <- resid(M1, type = "pearson")
N <- nrow(BF)
p <- length(fixef(M1)) + 1
sum(E1^2) / (N - p)
#Overdispersed!

#

#Fit the same model in glmmADMB
library(glmmADMB)
M2 <- glmmadmb(NEggs ~ NLowBranchesc + TreeHeightc + HerbCoverc + 
                    BuckthornAbundancec + BrambleAbundancec + 
                    SmallOakAbundancec + Distance2Edgec + NearestTallOakc +
                    fShelter +
                    (1 | fLandscapeID),
          family = "poisson",
          zeroInfl = FALSE,
          data = BF)
summary(M2)

#Check for overdispersion
E2 <- resid(M2, type = "pearson")
N <- nrow(BF)
p <- length(fixef(M2)) + 1 
sum(E2^2) / (N - p)

#It seems these results are different. Go for glmmadmb?


#Why is there overdispersion?
#Apply a detailed model validation


#Plot the Pearson residuals vs each covariate
MyVar <- c("NLowBranches",
          "TreeHeight", "HerbCover",
          "BuckthornAbundance", "BrambleAbundance", "SmallOakAbundance",   
          "Distance2Edge", "NearestTallOak", "X",  "Y" )

BF$E2 <- E2
Myxyplot(BF, MyVar, "E2", MyYlab = "Pearson residuals")

#Or the ggplot2 version of this graph:
MyVarx <- c("NLowBranches",
           "TreeHeight", "HerbCover",
           "BuckthornAbundance", "BrambleAbundance", "SmallOakAbundance",   
           "Distance2Edge", "NearestTallOak", "X",  "Y" )

MyVary <- "E2"
MyMultipanel.ggp2(BF, MyVarx, MyVary, 
                  ylab = "Pearson residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
#No scope for GAMs
##################

#But there seems to be 1 outlier!
par(mfrow = c(1,1), mar = c(5,5,3,3), cex.lab = 1.5)
plot(x = fitted(M2),
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline (h = 0, lty = 2)


#Use the identify command to figure out who the outlier is
#Either delete it (and rerun the analysis), or go on.


#I guess the large number of zeros is causing the overdispersion
#We can either go the zero inflated route, or we can apply a NB GLMM
#ZIP is not part of this course....hence NB GLMM it is.

#Aply a NB GLMM
M4 <- glmmadmb(NEggs ~ NLowBranchesc + TreeHeightc + HerbCoverc + 
                    BuckthornAbundancec + BrambleAbundancec + 
                    SmallOakAbundancec + Distance2Edgec + NearestTallOakc +
                    fShelter +
                    (1 | fLandscapeID),
          family = "nbinom",
          zeroInfl = FALSE,
          data = BF)
summary(M4)

#Check for overdispersion
E4 <- resid(M4, type = "pearson")
N <- nrow(BF)
p <- length(fixef(M4)) + 1 + 1
sum(E4^2) / (N - p)


#Lazy option: observation level random effect
#Not good if there are lots of zeros...which is the case!
BF$Eps <- 1 : nrow(BF)
BF$fEps <- factor(BF$Eps)
M6 <- glmmadmb(NEggs ~ NLowBranchesc + TreeHeightc + HerbCoverc + 
                    BuckthornAbundancec + BrambleAbundancec + 
                    SmallOakAbundancec + Distance2Edgec + NearestTallOakc +
                    fShelter +
                    (1 | fLandscapeID) + (1 | fEps),
          family = "poisson",
          zeroInfl = FALSE,
          data = BF)
summary(M6)


AIC(M1, M2, M4, M6)
   # df      AIC
# M1 12 744.3095   #Poisson GLMM in glmer
# M2 12 745.4060   #Poisson GLMM in glmmadmb
# M4 13 653.7160   #NB GLMM in glmmadmb
# M6 13 654.5500   #Poisson + observation level random intercept in glmmadmb



#Interpretation
summary(M4)
#What does it tell us??

#And is shelter significant?
drop1(M4, test = "Chi")
# Single term deletions

# Model:
# NEggs ~ NLowBranchesc + TreeHeightc + HerbCoverc + BuckthornAbundancec + 
    # BrambleAbundancec + SmallOakAbundancec + Distance2Edgec + 
    # NearestTallOakc + fShelter + (1 | fLandscapeID)
                    # Df    AIC    LRT  Pr(>Chi)    
# <none>                 653.72                     
# NLowBranchesc        1 661.40  9.688  0.001855 ** 
# TreeHeightc          1 653.67  1.954  0.162156    
# HerbCoverc           1 656.49  4.772  0.028926 *  
# BuckthornAbundancec  1 652.87  1.158  0.281881    
# BrambleAbundancec    1 652.09  0.374  0.540832    
# SmallOakAbundancec   1 652.74  1.028  0.310629    
# Distance2Edgec       1 680.60 28.886 7.677e-08 ***
# NearestTallOakc      1 652.09  0.374  0.540832    
# fShelter             2 655.19  5.470  0.064894 . 

#Model validation on M4.
E4 <- resid(M4, type = "pearson")
F4 <- fitted(M4)

par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = F4, 
     y = E4,
     xlab = "Fitted values (with re)",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)     


#How good is the model performing?
par(mfrow = c(1,1), cex.lab = 1.5, mar = c(5,5,2,2))
plot(x = BF$NEggs,
     y = F4,
     xlab = "Fitted values (with re)",
     ylab = "Number of eggs")
abline(h = 0, lty = 2)     

cor(BF$NEggs, F4)
source(file = "/Users/Highstat/applicat/HighlandStatistics/Courses/Data/R2glmms.R")
r.squared.lme(M4)  #Hmm..doesn't work for a NB GLMM via glmmadmb
#######################


#Dependency
MyCex <- abs(E4) / max(abs(E4)) 
MyCol <- sign(E4) + 2
xyplot(Y ~ X,
       data = BF,
       cex = 4 * sqrt(MyCex),
       col = MyCol)
# Homework: Make a variogram of the residuals!

#More dependency
MyVar <- c("NLowBranches",
          "TreeHeight", "HerbCover",
          "BuckthornAbundance", "BrambleAbundance", "SmallOakAbundance",   
          "Distance2Edge", "NearestTallOak", "X",  "Y" )

BF$E4 <- E4
Myxyplot(BF, MyVar, "E4", MyYlab = "Pearson residuals")

###########################################################


#Task:
#1. Write down the estimated model
#2. Make a picture shopwing the effects of the 
#   most important covariates 


