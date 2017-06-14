# this is the first code for modelling our data with lme4
# will incorporate RStan model as well


### INITIALIZE

library(lme4)
library(arm)
library(rstan)

# here set the working directory that points to the data folder
# e.g. the folder with annotated data saved as "myFA.Rdata"
# please just comment out the directories already here, so each
# user can uncomment the setwd for them as we move code around

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
setwd("D:/DataPlus2017/Data")


### LOAD DATA

# load fully annotated data (saved from LoadDataAndQC.R)
load("myFA.Rdata")

# extract CpG site 1 to start
firstSite <- head(FullAnnotation, 1)

# convert first site to table, using only H E K* W C (tumor+normal)
# 15 rows, with columns for patient, tumor/normal
# index in order A,B,N
indices <- c(11, 12, 13, 17, 18, 16, 20, 21, 19, 25, 26, 24, 33, 34, 32, 
             39, 31, 15)
patientLabel <- c(rep("H", 3),rep("E", 3),rep("K*", 3),rep("W", 3),rep("C", 3),
                  rep("J", 3))
tissueLabel <- rep(c("T", "T", "N"), 6)
sideLabel <- rep(c("A", "B", "NA"), 6)

firstData <- data.frame(beta=t(firstSite[indices]), patient = patientLabel, tissue=tissueLabel, side=sideLabel)


### FIT WITH LMER

fit1 <- lmer(formula = X1 ~ 1 + (1|patient/tissue), data=firstData)
summary(fit1)


### FIT WITH STAN

# insert fit with stan