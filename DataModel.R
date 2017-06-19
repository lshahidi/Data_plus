# this is the first code for modelling our data with lme4
# will incorporate RStan model as well


### INITIALIZE

# install.packages("coda")

library(lme4)
library(arm)
library(rstan)
library(coda)

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

firstData <- data.frame(beta=logit(t(firstSite[indices])), patient = patientLabel, tissue=tissueLabel, side=sideLabel)


# here we use all patient samples, excluding glands
indices <- c(9:13, 15:39,46,47,57,58,72:75)
patientLabel <- substr(colnames(firstSite[indices]),1,1)
patientLabel[10:12] <- "K*"
sideLabel <- substr(colnames(firstSite[indices]),2,2)
tissueLabel <- sideLabel
tissueLabel[tissueLabel %in% c("A","B")] <- "T"   #replace A and B with T
tumorIndicator <- 1*(tissueLabel=="T")

firstData <- data.frame(beta=logit(t(firstSite[indices])), patient = patientLabel, tissue=tissueLabel, side=sideLabel, tInd=tumorIndicator)


# here we use only patients with glands
indices <- c(9:13,32:34,57,58,60:71,76:79)
patientLabel <- substr(colnames(firstSite[indices]),1,1)
sideLabel <- substr(colnames(firstSite[indices]),2,2)
tissueLabel <- sideLabel
tissueLabel[tissueLabel %in% c("A","B")] <- "T"   #replace A and B with T
tumorIndicator <- 1*(tissueLabel=="T")
glandIndicator <- 1*(indices %in% c(60:71,76:79))

firstData <- data.frame(beta = logit(t(firstSite[indices])), patient = patientLabel, tissue=tissueLabel, side=sideLabel, tInd=tumorIndicator, gInd=glandIndicator)




### FIT WITH LMER

fit1 <- lmer(X1 ~ tInd + (tInd|patient), firstData)

fitGland <- lmer(X1 ~ tInd + (tInd|patient) + (0+gInd|patient:side), firstData)

patLabs <- c("C", "D","E","F","H","J","K","K*","M","O","P","S","T","U","W","X")
# barplot(fixef(fit1), main="Fixed Effects", xlab="Effect", ylab="Intercept Estimate", names.arg=c("mu","betaT"))
# barplot(ranef(fit1)$patient$'(Intercept)', main="Patient Random Effect", xlab="Patient", ylab=expression(paste("Intercept Estimate (",b[i],")")), cex.names = 0.7, las=2)
# barplot(ranef(fit1)$patient$tInd, main="Patient,Tumor Slope Random Effect", xlab="Patient", ylab=expression(paste("Slope Estimate (",b[iT],")")), cex.names = 0.7, las=2)

barplot(fixef(fit1), main="Fixed Effects", xlab="Effect", ylab="Intercept Estimate", names.arg=c("mu","betaT"), ylim=c(-0.2,0.2))
barplot(fixef(fit2), main="Fixed Effects", xlab="Effect", ylab="Intercept Estimate", names.arg=c("mu","betaT"), ylim=c(-0.2,0.2))
barplot(ranef(fit3)$patient$'(Intercept)', main="Patient Random Effect", xlab="Patient", ylab=expression(paste("Intercept Estimate (",b[i],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-0.2,0.2))
barplot(ranef(fit2)$patient$'(Intercept)', main="Patient Random Effect", xlab="Patient", ylab=expression(paste("Intercept Estimate (",b[i],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-0.2,0.2))
barplot(ranef(fit3)$patient$tInd, main="Patient,Tumor Slope Random Effect", xlab="Patient", ylab=expression(paste("Slope Estimate (",b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-0.2,0.2))
barplot(ranef(fit2)$patient$tInd, main="Patient,Tumor Slope Random Effect", xlab="Patient", ylab=expression(paste("Slope Estimate (",b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-0.2,0.2))

barplot(ranef(fit1)$patient$'(Intercept)'+ranef(fit1)$patient$tInd, main="Both Random Effects", xlab="Patient", ylab=expression(paste("Both Estimates (",b[i]+b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-0.2,0.2))
barplot(ranef(fit2)$patient$'(Intercept)'+ranef(fit2)$patient$tInd, main="Both Random Effects", xlab="Patient", ylab=expression(paste("Both Estimates (",b[i]+b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-0.2,0.2))

#barplot(fixef(fit1), main="Patient Fixed Effects", xlab="Patient", ylab=expression(paste("Intercept Estimate (",beta[k],")")), names.arg=c("C", "D","E","F","H","J","K","K*","M","O","P","S","T","U","W","X"), cex.names=0.7)
#barplot(fixef(fit1)[-17]+fixef(fit1)[1], main="Patient Fixed Effects", xlab="Patient", ylab=expression(paste("Intercept Estimate (",beta[k],")")), names.arg=c("C","D","E","F","H","J","K","K*","M","O","P","S","T","U","W","X"), cex.names=0.65)




### FIT WITH STAN

# Create tumor indicator

for (i in 1:dim(firstData)[1]) {
  if (firstData$tissue[i] == "T") {
    firstData$tInd[i] = 1
  } else {
    firstData$tInd[i] = 0
  }
}

# Logit transform on beta values

#firstData$y <- log(firstData$X1/(1-firstData$X1))
firstData$y <- firstData$X1

# Stan data

stanDat <- list(pID = as.integer(factor(firstData$patient)),
                tInd = firstData$tInd,
                N = nrow(firstData),
                P = nlevels(firstData$patient),
                y = firstData$y)


stanFit <- stan(file="model.stan",data=stanDat)

print(stanFit, probs = c(0.025, 0.5, 0.975))

summary(stanFit)

est <- as.matrix(stanFit,pars=c("mu","betaT","b_pat","bT_pat"))
# Estimates for first three samples
head(est[,c(1:5,19:21)])


# Diagnostic Graphs
mcmcCoda <- mcmc.list(lapply( 1:ncol(stanFit) ,
                              function(x) { mcmc(as.array(stanFit)[,x,])} ))

plot(stanFit, pars=c("b_pat"), main="Random Effect of Patients")
plot(stanFit, pars=c("bT_pat"), main="Random Effect of Tumor Tissue Varying Between Patients")



