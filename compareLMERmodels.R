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

# extract variance at 1000 sites
nSites<-50
sigmaLMER1 <- data.frame(sigmaT=numeric(nSites),sigmaP=numeric(nSites),sigmaE=numeric(nSites))
sigmaLMER2 <- data.frame(sigmaT=numeric(nSites),sigmaP=numeric(nSites),sigmaE=numeric(nSites))

ptm <- proc.time()
for (i in 1:nSites){
  print(paste("site ",i))

  DataI <- site(i)
  names(DataI)[1] <- "beta"

  fit1 <- lmer(beta ~ tInd + (tInd|patient), DataI)

  sigmaLMER1$sigmaP[i] <- as.data.frame(VarCorr(fit1))$sdcor[1]
  sigmaLMER1$sigmaT[i] <- as.data.frame(VarCorr(fit1))$sdcor[2]
  sigmaLMER1$sigmaE[i] <- as.data.frame(VarCorr(fit1))$sdcor[4]
  sigmaLMER1$betaT[i] <- fixef(fit1)[2]
  sigmaLMER1$avT[i] <- mean(ranef(fit2)$patient$tInd)
  sigmaLMER1$mu[i] <- fixef(fit1)[1]
  
  fit2 <- lmer(beta ~ (tInd|patient), DataI)
  
  sigmaLMER2$sigmaP[i] <- as.data.frame(VarCorr(fit2))$sdcor[1]
  sigmaLMER2$sigmaT[i] <- as.data.frame(VarCorr(fit2))$sdcor[2]
  sigmaLMER2$sigmaE[i] <- as.data.frame(VarCorr(fit2))$sdcor[4]
  sigmaLMER2$avT[i] <- mean(ranef(fit2)$patient$tInd)
  sigmaLMER2$mu[i] <- fixef(fit2)[1]
}
proc.time() - ptm