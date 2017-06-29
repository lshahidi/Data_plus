# compare lmer models with and without fixed effects

### INITIALIZE

# libraries
library(lme4)
library(arm)
library(gtools)
library(reshape2)
library(ggplot2)

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
setwd("D:/DataPlus2017/Data")

# load data
load("myFA.Rdata")

site <- function (site_no) {
  # extract CpG site xx to start
  temp <- FullAnnotation[site_no,]
  
  # here we use all patient samples, excluding glands
  indices <- c(9:13, 15:39,46,47,57,58,72:75)
  patientLabel <- substr(colnames(temp[indices]),1,1)
  patientLabel[10:12] <- "K*"
  sideLabel <- substr(colnames(temp[indices]),2,2)
  tissueLabel <- sideLabel
  tissueLabel[tissueLabel %in% c("A","B")] <- "T"   #replace A and B with T
  tumorIndicator <- 1*(tissueLabel=="T")
  
  work <- data.frame(beta=logit(t(temp[indices])), patient = patientLabel, tissue=tissueLabel, side=sideLabel, tInd=tumorIndicator)
  
  return(work)
}

siteNorm <- function (site_no) {
  # extract CpG site xx to start
  temp <- FullAnnotation[site_no,]
  
  # here we use all patient samples, excluding glands
  indices <- c(11:13,15:21,24:26,31:34,39)
  patientLabel <- substr(colnames(temp[indices]),1,1)
  patientLabel[8:10] <- "K*"
  sideLabel <- substr(colnames(temp[indices]),2,2)
  tissueLabel <- sideLabel
  tissueLabel[tissueLabel %in% c("A","B")] <- "T"   #replace A and B with T
  tumorIndicator <- 1*(tissueLabel=="T")
  
  work <- data.frame(beta=logit(t(temp[indices])), patient = patientLabel, tissue=tissueLabel, side=sideLabel, tInd=tumorIndicator)
  
  return(work)
}

# extract variance at 1000 sites
nSites<-500
sigmaLMER1 <- data.frame(sigmaT=numeric(nSites),sigmaP=numeric(nSites),sigmaE=numeric(nSites))
sigmaLMER2 <- data.frame(sigmaT=numeric(nSites),sigmaP=numeric(nSites),sigmaE=numeric(nSites))

ptm <- proc.time()
for (i in 1:nSites){
  print(paste("site ",i))

  DataI <- site(i)
  names(DataI)[1] <- "logitbeta"

  fit1 <- lmer(logitbeta ~ tInd + (tInd||patient), DataI)

  sigmaLMER1$sigmaP[i] <- as.data.frame(VarCorr(fit1))$sdcor[1]
  sigmaLMER1$sigmaT[i] <- as.data.frame(VarCorr(fit1))$sdcor[2]
  sigmaLMER1$sigmaE[i] <- as.data.frame(VarCorr(fit1))$sdcor[3]
  sigmaLMER1$betaT[i] <- fixef(fit1)[2]
  sigmaLMER1$avT[i] <- mean(ranef(fit2)$patient$tInd)
  sigmaLMER1$mu[i] <- fixef(fit1)[1]
  
  DataEst = data.frame(patient=rep(c("C","D","E","F","H","J","K","K*","M","O","P","S","T","U","W","X"),2), tInd=c(rep(0,16),rep(1,16)))
  mu <- fixef(fit1)[1]
  betaT <-fixef(fit1)[2]
  bi <- ranef(fit1)$patient[DataEst$patient,1]
  bTi <- ranef(fit1)$patient[DataEst$patient,2]
  sigmaP <- as.data.frame(VarCorr(fit1))$sdcor[1]
  sigmaT <- as.data.frame(VarCorr(fit1))$sdcor[2]
  sigmaE <- as.data.frame(VarCorr(fit1))$sdcor[3]
  DataEst$estBeta <- mu + betaT*DataEst$tInd + bi + bTi*DataEst$tInd
  print(ggplot() + geom_point(data=DataI, aes(x = tInd, y=logitbeta, color=patient, group=patient)) + geom_line(data=DataEst, aes(x = tInd, y=estBeta, color=patient, group=patient)) + labs(title=paste("Site ",as.integer(i))) + annotate("text", x=c(0.25, 0.75, 0.5), y = yHeights, label=c(paste("sigmaP: ",sigmaP),paste("sigmaT: ",sigmaT),paste("sigmaE: ",sigmaE))))
  readline(prompt="Press [enter] to continue")
}
proc.time() - ptm


# using uncorrelating model to fit
siteNum <- 5
yHeights = c(-0.2, -0.2, -0.5)
DataI <- site(siteNum)
names(DataI)[1] <- "logitbeta"
fit <- lmer(logitbeta ~ tInd + (tInd||patient), DataI)
# estimate beta for all 16 patients at tInd=[0,1]
DataEst = data.frame(patient=rep(c("C","D","E","F","H","J","K","K*","M","O","P","S","T","U","W","X"),2), tInd=c(rep(0,16),rep(1,16)))
mu <- fixef(fit)[1]
betaT <-fixef(fit)[2]
bi <- ranef(fit)$patient[DataEst$patient,1]
bTi <- ranef(fit)$patient[DataEst$patient,2]
sigmaP <- as.data.frame(VarCorr(fit))$sdcor[1]
sigmaT <- as.data.frame(VarCorr(fit))$sdcor[2]
sigmaE <- as.data.frame(VarCorr(fit))$sdcor[3]
DataEst$estBeta <- mu + betaT*DataEst$tInd + bi + bTi*DataEst$tInd
ggplot() + geom_point(data=DataI, aes(x = tInd, y=logitbeta, color=patient, group=patient)) + geom_line(data=DataEst, aes(x = tInd, y=estBeta, color=patient, group=patient)) + labs(title=paste("Site ",as.integer(siteNum))) + annotate("text", x=c(0.25, 0.75, 0.5), y = yHeights, label=c(paste("sigmaP: ",sigmaP),paste("sigmaT: ",sigmaT),paste("sigmaE: ",sigmaE)))


# same but for correlating model
yHeights = c(-0.2, -0.2, -0.9)
DataI <- site(siteNum)
names(DataI)[1] <- "logitbeta"
fit <- lmer(logitbeta ~ tInd + (tInd|patient), DataI)
# estimate beta for all 16 patients at tInd=[0,1]
DataEst = data.frame(patient=rep(c("C","D","E","F","H","J","K","K*","M","O","P","S","T","U","W","X"),2), tInd=c(rep(0,16),rep(1,16)))
mu <- fixef(fit)[1]
betaT <-fixef(fit)[2]
bi <- ranef(fit)$patient[DataEst$patient,1]
bTi <- ranef(fit)$patient[DataEst$patient,2]
sigmaP <- as.data.frame(VarCorr(fit))$sdcor[1]
sigmaT <- as.data.frame(VarCorr(fit))$sdcor[2]
sigmaE <- as.data.frame(VarCorr(fit))$sdcor[4]
DataEst$estBeta <- mu + betaT*DataEst$tInd + bi + bTi*DataEst$tInd
ggplot() + geom_point(data=DataI, aes(x = tInd, y=logitbeta, color=patient, group=patient)) + geom_line(data=DataEst, aes(x = tInd, y=estBeta, color=patient, group=patient)) + labs(title=paste("Site ",as.integer(siteNum))) + annotate("text", x=c(0.25, 0.75, 0.5), y = yHeights, label=c(paste("sigmaP: ",sigmaP),paste("sigmaT: ",sigmaT),paste("sigmaE: ",sigmaE)))
