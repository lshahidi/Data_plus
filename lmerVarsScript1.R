# script for extracting variance from lmer model 
# part 1: on sites 1-425K

### INITIALIZE

# libraries
library(lme4)
library(arm)
library(gtools)
library(reshape2)
library(ggplot2)

# load data
load("myFA.Rdata")

# Function used to read in data from each site

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

### FIT WITH LMER

## extract variance at all sites
nSites<-dim(FullAnnotation)[1]
sigmaLMER <- data.frame(sigmaT=numeric(nSites),sigmaP=numeric(nSites),sigmaE=numeric(nSites))

ptm <- proc.time()
for (i in 1:nSites){
  print(paste("site ",i))
  
  DataI <- site(i)
  names(DataI)[1] <- "beta"
  
  fit <- lmer(beta ~ tInd + (tInd|patient), DataI)
  
  sigmaLMER$sigmaT[i] <- as.data.frame(VarCorr(fit))$sdcor[1]
  sigmaLMER$sigmaP[i] <- as.data.frame(VarCorr(fit))$sdcor[2]
  sigmaLMER$sigmaE[i] <- as.data.frame(VarCorr(fit))$sdcor[4]
}
proc.time() - ptm

# save data
save(sigmaLMER,file="myLMERVars.RData")