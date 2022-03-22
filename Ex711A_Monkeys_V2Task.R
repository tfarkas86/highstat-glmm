#    Highland Statistics Ltd.
#    www.highstat.com
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

#This exercise is based on:
# Evidence for varying social strategies across the day in chacma baboons
# Claudia Sick, Alecia J. Carter, Harry H. Marshall, Leslie A. Knapp, 
# Torben Dabelsteen and Guy Cowlishaw
# Biol. Lett. 2014 10, 20140249, published 9 July 2014
# Data where taken from the online supplemental material:
# http://rsbl.royalsocietypublishing.org/content/suppl/2014/07/08/rsbl.2014.0249.DC1.ht ml


#Description of experiment/data:
 #1. Data were collected from two troops of chacma baboons
 #2. Each troop was followed daily from dawn to dusk (ca 06.00–18.00 h). 
 #3. 1 h long focal follows on all identifiable individuals 
 #4. Recording all grooming and dominance interactions
 #5. Unit of analysis:  grooming session
 #6. Multiple observations from the same focal session
 #7. Multiple observations on the same animal.
   


#What will you learn in this exercise?
 # 1A. Reckognising where the correlation is (2 way nested and crossed data)
 # 1B. Some simple ggplot2 coding
 # 1C. Binomial glmer coding
 # 1D. R^2 for a GLMM


###########################################################################
#Set the working directory
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/FilesMCMC_GLMM/CDMCMCGLMGLMMCourse/Data")	
Monkeys <- read.table(file = "MonkeysV2.txt", header = TRUE)
names(Monkeys)


#########################################################################################
# In the Gaussian mixed modelling exercises we analysed the 
# Monkey data. The question was: 
# Subordinates should prefer to groom more dominant animals 
# earlier in the day
# And we used these data to answer the question:

Monkeys.sub <- Monkeys[Monkeys$SubordinateGrooms == "yes",]
MS <- Monkeys.sub  #Is shorter

# Now we focus on a different question for this data set.
# There is also a variable GroomSymmetry:
#    -binary term, where 1 indicated that the subordinate 
#     groomed more than 50% of the session’s full duration.

#(I am not sure whether GroomSymmetry is the correct word here).


# Recall the setup of the experiment
# 1. Two monkey groups (large vs small). Coded as GroupSize
# 2. There are multiple focal hours per day (coded as "FocalHour")
# 3. In a focal hour we look for a groomer (coded as FocalGroomer)
# 4. We also write down the id of the receiver monkey  (coded as Receiver)
# 5. Record the time since sunrise, and the rank difference between 
#    focal groomer and receiver

#Aim: GroomSymmetry = function(Time, Relatedness, GroupSize)

#Interactions are expected:
    #Time and relatedness  (used in paper)
    #Based on our common sense: Relatedness and groupsize  ?
############################################################################



############################################################################
#Understanding the dependency in the data
#Where is the correlation?
#Our first idea was:
#All observations made in the same hour
#All observations made on the same focal groomer
#All observations made on the same receiver

#This would mean:
#(1 | FocalHour) + (1 | FocalGroomer) + (1 | Receiver) 



###########################################################################
#Response variable: GroomSymmetry  (0/1)
#Covariates: Time           (time since sunrise)
#            Relatedness    (index specifying relatedness)
#            GroupSize      (small vs large)  
#
#            FocalGroomer     (monkey id who is doing the grooming)
#            FocalHour        (focal hour)
#            Receiver         (grooming partner identity)
#
#######################################################################################




#######################################################################################
#Housekeeping
library(lattice)
source(file = "/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV9.R")
library(ggplot2)
library(lme4)
#################################################################

#Task:
# 1. Apply data exploration (condensed)
# 2. Model GroomSymmetry as a function of the covariates:
#         Time * Relatedness + Relatedness * GroupSize +
#         (1 | FocalGroomer) + (1 | FocalHour) + (1 | Receiver)
#    Use frequentist tools
# 3. Apply model selection
# 4. Sketch the fitted model, and write down the equation
# 5. Try getting an R2.

#################################################################

