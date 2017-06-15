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

# convert first site to table, using only H E K* W C J (tumor+normal)
# 18 rows, with columns for patient, tumor/normal
# index in order A,B,N
indices <- c(11, 12, 13, 17, 18, 16, 20, 21, 19, 25, 26, 24, 33, 34, 32, 
             39, 31, 15)
patientLabel <- c(rep("H", 3),rep("E", 3),rep("K*", 3),rep("W", 3),rep("C", 3),
                  rep("J", 3))
tissueLabel <- rep(c("T", "T", "N"), 6)
sideLabel <- rep(c("A", "B", ""), 6)

# here we use all samples
indices <- c(9:13, 15:39,46,47,57,58,72:75)
patientLabel <- substr(colnames(firstSite[indices]),1,1)
patientLabel[10:12] <- "K*"
sideLabel <- substr(colnames(firstSite[indices]),2,2)
tissueLabel <- sideLabel
tissueLabel[tissueLabel %in% c("A","B")] <- "T"   #replace A and B with T

firstData <- data.frame(beta=t(firstSite[indices]), patient = patientLabel, tissue=tissueLabel, side=sideLabel)


### FIT WITH LMER

fit1 <- lmer(X1 ~ patient + (1|patient:tissue), firstData)
barplot(fixef(fit1)+fixef(fit1)[1], main="Patient Fixed Effects", xlab="Patient", ylab=expression(paste("Intercept Estimate (",beta[k],")")), names.arg=c("C","D","E","F","H","J","K","K*","M","O","P","S","T","U","W","X"), cex.names=0.65)
barplot(t(as.vector(ranef(fit1)$'patient:tissue')), main="Patient:Tissue Random Effects", xlab="Patient:Tissue", ylab=expression(paste("Intercept Estimate (",epsilon[j]^alpha,")")), cex.names = 0.7, las=2)

barplot(fixef(fit1), main="Patient Fixed Effects", xlab="Patient", ylab=expression(paste("Intercept Estimate (",beta[k],")")), names.arg=c("C", "D","E","F","H","J","K","K*","M","O","P","S","T","U","W","X"), cex.names=0.7)

### FIT WITH STAN

# insert fit with stan