# script for estimating sigma priors from TCGA patient data
# loads and combines patient data (with annotation)
# then algebraically calculates prior distribution estimates for each site
# to ultimately produce a full prior distribution then fit with a parameterized gamma
# currently test version

library(gtools)
library(ggplot2)

### INITIALIZE ###

log.file <- "logTCGA.txt"
writeLines(c(""), log.file)

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/TCGA data")

# Yanlin's working directory
#setwd("D:/DataPlus2017/Data")

# Load and concatenate all available data
# need to: search all folders for data, load each one just the beta and label
fileList <- list.files(pattern = "_hg38.txt$", recursive = TRUE)

# intialize data frame
ptmTot <- proc.time()
print("Initializing data frame with annotations")
temp <- read.table(fileList[1],header=TRUE,sep="\t",fill=TRUE)
concatData <- as.data.frame.matrix(temp[,-2])

print("Now adding patient data")
for(f in fileList) {
  ptm <- proc.time()
  print(paste("File: ",which(f==fileList),"/",length(fileList)))
  dataLabel <- substr(f,regexpr("TCGA",f),regexpr("TCGA",f)+15)
  print(paste("Label: ", dataLabel))
  temp <- read.table(f,header=TRUE,sep="\t",fill=TRUE)
  concatData[dataLabel] <- temp[,2]
  print(proc.time() - ptm)
}
rm(temp)
gc()

print("Total time to load data:")
print(proc.time()-ptmTot)

save(concatData, "testdata3pats.Rdata")
load("testdata3pats.Rdata")

# create labels
# start by making vectors for labels of patient and normal vs tumor
dataLabel <- colnames(concatData)[-1:-10]
patLabel <- substr(dataLabel,9,12)
tumorLabel <- 1-as.numeric(substr(dataLabel,14,14))
sideLabel <- substr(dataLabel,16,16)
sideLabel[grep("C",sideLabel)] <- "B" # simply re-labels any C tumor with B (should be only one case)
sideLabel[!tumorLabel] <- "N"

# find sites with NA
# need to check that all samples are not NA? or skip with any NA? or maybe at least half
# currently going to remove any site with NA
NAsites <- c()
for(i in dataLabels){
  NAsites <- union(NAsites,which(is.na(concatData[,i])))
}
print(paste("# sites with at least one NA: ",length(NAsites)))

 # randomly select 5 sites in 1-485577 (excluding sites with any NA samples)
nsites <- 100
siteInds <- sample(setdiff(1:485577,NAsites),nsites)


### FXNS ###
# Extract single site data
site <- function (site_no) {
  # extract CpG site xx to start, exclude first 10 columns
  temp <- concatData[site_no,-1:-10]
  
  # already have dataLabel, patLabel, tumorLabel, sideLabel defined in initialization
  tissueLabel <- sideLabel
  tissueLabel[tissueLabel %in% c("A","B")] <- "T"   #replace A and B with T
  tumorIndicator <- 1*(tissueLabel=="T")
  
  work <- data.frame(beta=logit(t(temp)), patient = patLabel, tissue=tissueLabel, side=sideLabel, tInd=tumorIndicator)
  colnames(work)[1] <- "beta"
  return(work)
}

# Plot raw data and fits
fitplot <- function(data, line) {
  # plot beta and fits at this site for all patients, using mean estimates for b and c coefficients
  p1 <- ggplot() + geom_point(data = data, aes(x=tInd, y=beta, colour=patient))
  pats <- as.integer(factor(data$patient))
  npats <- max(pats)
  mu <- line[1]
  betaT <- line[2]
  linedf <- data.frame(tInd = c(0,1), est = c(mu, (mu+betaT)))
  p1 <- p1 + geom_line(data=linedf, aes(x=tInd, y=est))
  
  return(p1)
}

# function to calculate sigmaP, uses only normal samples
# mu is mean, b is diff from mu, sigmaP is sample variance
calcSigmaP <- function (siteBeta) {
  mu <- mean(siteBeta[!tumorLabel])
  b <- siteBeta[!tumorLabel] - mu
  print(paste("length b: ", length(b)))
  sigmaP <- var(b)
  return(c(mu,sigmaP))
}

# function to calculate sigmaPT, uses normal and tumor samples
# first take tumor mean (if multiple), then normalize by subtracting patient normal
# then mean is betaT, c is diff from betaT, sigmaPT is sample variance
calcSigmaPT <- function (siteBeta) {
  normBeta <- c()
  for(pat in unique(patLabel)){
    if(length(intersect(which(pat==patLabel),which(!tumorLabel)))){ # only if pat has normal
      normIndex <- intersect(which(pat==patLabel),which(!tumorLabel)) # normal sample index
      tumIndex <- intersect(which(pat==patLabel),which(!!tumorLabel)) # tumor sample indices
      if(length(tumIndex)-1){ tum<-mean(siteBeta[tumIndex]) }else{ tum<-siteBeta[tumIndex] }
      normBeta <- c(normBeta, tum-siteBeta[normIndex])
    }
  }
  betaT <- mean(normBeta)
  c <- normBeta - mu
  print(paste("length c: ", length(c)))
  sigmaPT <- var(c)
  return(c(betaT,sigmaPT))
}

# function to calculate sigmaT, uses patients with multiple tumor samples
# normalize by subtracting patient tumor mean
# then these values are d, sigmaT is sample variance
calcSigmaT <- function (siteBeta) {
  normBeta <- c()
  for(pat in unique(patLabel)){
    if(length(intersect(which(pat==patLabel),which(!!tumorLabel)))-1){ # only if pat has multiple tumor
      tumIndex <- intersect(which(pat==patLabel),which(!!tumorLabel)) # tumor sample indices
      mean <- siteBeta[tumIndex]
      normBeta <- c(normBeta, siteBeta[tumIndex] - mean)
    }
  }
  d <- normBeta
  print(paste("length d: ", length(d)))
  sigmaT <- var(d)
  return(sigmaT)
}


### CODE ###
# intialize numeric vectors for each site estimate
sigmaP_ests <- numeric(dim(concatData)[1])
sigmaPT_ests <- sigmaP_ests
sigmaT_ests <- sigmaP_ests
mu_ests <- sigmaP_ests
betaT_ests <- sigmaP_ests

ptm<-proc.time()
for(i in siteInds){
  print(paste("Analyzing site ",i,"."))
  siteData <- site(i)
  
  temp <- calcSigmaP(siteData$beta)
  mu_ests[i] <- temp[1]
  sigmaP_ests[i] <- temp[2]
  temp <- calcSigmaPT(siteData$beta)
  betaT_ests[i] <- temp[1]
  sigmaPT_ests[i] <- temp[2]
  sigmaT_ests[i] <- calcSigmaT(siteData$beta)
  print(paste("Site ",i," completed."))
  
  # plot raw data with mu and beta
  #p1 <- fitplot(site(i),c(mu_ests[i],betaT_ests[i]))
  #print(p1)
  
  #invisible(readline(prompt="Press [enter] to continue"))
}
print(proc.time()-ptm)
