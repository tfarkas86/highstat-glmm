
#    Highland Statistics Ltd.
#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


######################################################################
# Mallard exercise
# These data are taken from:
# Transfer of Maternal Antibodies against Avian Influenza# Virus in Mallards (Anas platyrhynchos)# Jacintha G. B. van Dijk, A. Christa Mateman, Marcel Klaassen
# PLOS ONE November 2014 | Volume 9 | Issue 11 | e112595

#Set the working directory and import the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/NewDataSets/lme_mallards/")
Mallards <- read.table("Mallards.txt", 
                   header = TRUE,
                   dec = ".")
names(Mallards)
str(Mallards)

# From the abstract
# Maternal antibodies protect chicks from infection 
# with pathogens early in life and may impact pathogen 
# dynamics due to the alteration of the proportion of 
# susceptible individuals in a population. We investigated 
# the transfer of maternal antibodies against avian influenza 
# virus (AIV) in a key AIV host species, the mallard 
# (Anas platyrhynchos). Combining observations in both the 
# field and in mallards kept in captivity, we connected 
# maternal AIV antibody concentrations in eggs to
#   (i) female body condition, 
#  (ii) female AIV antibody concentration, 
# (iii) egg laying order,    #<--- but only in the capticity study
#  (iv) egg size and 
#   (v) embryo sex.

#Methods: 
# 1. FIELD STUDY DATA
# 2. Summer 2010: 67 free-living female mallards 
#    were caught from their nest with a sweep net
# 3. To index body size: we measured tarsus length, 
#    head+bill length and wing length. Plus body mass
# 4. Blood samples to detect antibodies to AIV.
# 5. Per clutch: two randomly chosen eggs were collected 
#    to assess maternal AIV antibody concentration in 
#    egg yolk. (=RESPONSE VARIABLE)
#6.  Quantified size/gender of embryo 

#Summary:
#  64 female nesting ducks. 2 eggs per duck/clutch
#  115 observations (there is always some sampling trouble)
###################################

# Task: Model Antibody concentrations in egg yolk as a 
#       function of:
#          -Body size variable (use one if collinear)
#          -Antibodies in mother duck
#          -Mass of mother duck
#          -EggVolume 
#          -EmbryoSize 
#          -EmbryoSex
#      But we have multiple obervations per clutch!
########################################################




########################################################
#Coding:
Mallards$abYolk  <- -Mallards$ODvalueeggyolk  #Response
Mallards$abMom   <- -Mallards$ODvaluefemale   #Covariate
Mallards$fClutch <- factor(Mallards$Clutch)   #Random effect

#The minus signs for antibodies in yolk and in the female
#was used in the paper. It makes the interpretation easier.
########################################################



###########################################################
#Load packages and library files
library(lattice)  #Fancy multi-panel graphs
library(nlme)     #Mixed modelling
library(mgcv)     #Smoothing...needed by HighstatLib
library(ggplot2)  #Fancy multi-panel graphs

#source the file: HighstatLibV7.R
source("/Users/Highstat/MaggieWindows/applicat/HighlandStatistics/webHighstatNew2_2010/CourseMCMCGLMM/RCode/HighstatLibV7.R")
##################################################################


#Task: Carry out a data exploration, frequentist analysis and Bayesian analysis 


