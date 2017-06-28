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


DataI <- site(1)
names(DataI)[1] <- "logitbeta"
fit <- lmer(logitbeta ~ tInd + (tInd|patient), DataI)
#fit <- lmer(logitbeta ~ tInd + (1|patient) + (0+tInd|patient), DataI)
mu <- fixef(fit)[1]
betaT <-fixef(fit)[2]
bi <- ranef(fit)$patient[DataI$patient,1]
bTi <- ranef(fit)$patient[DataI$patient,2]
sigmaP <- as.data.frame(VarCorr(fit))$sdcor[1]
sigmaT <- as.data.frame(VarCorr(fit))$sdcor[2]
sigmaE <- as.data.frame(VarCorr(fit))$sdcor[4]
DataI$estBeta <- mu + betaT*DataI$tInd + bi + bTi*DataI$tInd
ggplot() + geom_point(data=DataI, aes(x = tInd, y=logitbeta, color=patient, group=patient)) + geom_line(data=DataI, aes(x = tInd, y=estBeta, color=patient, group=patient)) + labs(title="Site 1") + annotate("text", x=c(0.25, 0.75, 0.5), y = c(-0.5, -0.5, -0.7), label=c(paste("sigmaP: ",sigmaP),paste("sigmaT: ",sigmaT),paste("sigmaE: ",sigmaE)))
