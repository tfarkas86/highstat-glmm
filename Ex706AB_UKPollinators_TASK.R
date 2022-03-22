#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#############################################################################
#UK pollinator data
#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/")

#Abundance data
PolAB <- read.table("UKPollinatorsV2_TotalAbundance.txt", 
                    header = TRUE,
                    dec = ".")
#Richness data
PolRI <- read.table("UKPollinatorsV2_Richness.txt", 
                    header = TRUE,
                    dec = ".")
#Covariates
Cov <- read.table("UKPollinatorsV2_Covariates.txt", 
                  header = TRUE,
                  dec = ".")
str(PolAB)
str(PolRI)
str(Cov)                  
#############################################################################


# Data from:
# Where is the UK’s pollinator biodiversity?# The importance of urban areas for flowervisiting# insects.
# Baldock et al. (2015). 
# Proc. R. Soc. B 282: 20142849.# http://dx.doi.org/10.1098/rspb.2014.2849


# From the abstract:
# Insect pollinators provide a crucial ecosystem service, 
# but are under threat. Urban areas could be important for 
# pollinators, though their value relative to other habitats 
# is poorly known. We compared pollinator communities# using quantified flower-visitation networks in 36 sites 
# (each 1 km2) in three landscapes: urban, farmland and 
# nature reserves. ....


# In short:
# 1. The 36 sites were located in and around 12 large 
#    UK urban centres with populations over 150 000. 
# 2. Cities were blocked into four regional groups of 
#    three (farmland, urban, nature reservate ). 
#		-Urban sites were located within the city boundary
#       -Farmland and nature reserve sites within 10 km of 
#        the city boundary.
# 3. Each of the 36 sites was sampled four times between 
#    30 May and 19 September 2011 at approximately monthly 
#    intervals (using transects).
# 4. Flowers were sampled at 10 m intervals along each transect.
# 5. Flower-visitor interactions were quantified by walking along#    each transect and collecting all insects on flowers up to 
#    1 m either side of the transect line to a height of 2 m.


# Response variable: Species richness and visitor abundance
# Covariates: 
#       - Landscape type (urban,farmland, nature reserve)
#       - Sampling month (June, July, August, September)
#       - Floral abundance 
#       - Proportion of woodland habitat
# Multiple observations per city!

# In the paper:
# Analyses were carried out for (i) the whole dataset; 
# (ii) separately for the two dominant insect orders, 
# Diptera and Hymenoptera; (iii) for the key pollinator# taxa of hoverflies (Diptera: Syrphidae) and bees (Apoidea:# comprising bumblebees, honeybees and solitary bees); and# (iv) separately for bumblebees, honeybees and solitary bees.

# In this exercise we will focus on the abundance of Diptera




######################################################################
#TASK: Model Diptera abundance as a function of Type, Month,
#      Floral abundance and proportion of woodland habitat.

#PRIME VARIABLE OF INTEREST: Type

#Use frequentist and Bayesian techniques

###################################################################
#Load packages and library files

library(lattice)
library(lme4) 
library(ggplot2) 
library(plyr)
source("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV9.R")


###############################################
#Housekeeping: Merge covariate with total abundance data

names(Cov)
PolAB$Woodland <- rep(Cov$Woodland_habitat, each = 4)
PolAB$Floral   <- rep(Cov$Floral_Unit_Abundance, each = 4)

#Because we may want to focus on different species,
#we put the selected data into a new data frame.
PolData <- data.frame(
            Y          = PolAB$Diptera,  #Species<--
            Type       = PolAB$Type,
            Month      = PolAB$Month,
            City       = PolAB$City,
            Woodland   = PolAB$Woodland,
            Floral     = PolAB$Floral,
            Woodland.c = MyStd(PolAB$Woodland),
            Floral.c   = MyStd(PolAB$Floral))

# Now you can easily change the species and rerun all 
# the code below without the need to change much.




######################################################################
#Data exploration
