# load StanCRun data
# produce full size data frames for mu,betaT,sigmas (P,PT,T,E)
# read in individual save files containing 1000 site runs
# store 10000 site chunks in corresponding positions in full data frame
# repeat for each completed run

library(ggplot2)
library(reshape2)

### LOADING FOR STAN MODEL RESULTS

# Set working directory to normal data wd then folder ResultsC
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data/ResultsC")

betaT_Cfull <-
  data.frame(
    mean = numeric(866836),
    p2.5 = numeric(866836),
    p25 = numeric(866836),
    p50 = numeric(866836),
    p75 = numeric(866836),
    p97.5 = numeric(866836)
  )

mu_Cfull <- betaT_Cfull
sigmaE_Cfull <- betaT_Cfull
sigmaP_Cfull <- betaT_Cfull
sigmaPT_Cfull <- betaT_Cfull
sigmaT_Cfull <- betaT_Cfull

nsites <- 10000
filesPresent <- (1:87)
fileNames <- dir(pattern = "StanCParResults_*")
allInds <- numeric()
for (i in filesPresent) {
  print(paste("Chunk:",i))
  load(fileNames[fileNames == paste("StanCParResults_",i,".Rdata", sep="")])
  inds <- (1:nsites) + (i - 1) * nsites
  if (i==87) {
    inds <- inds[1:6836]
  }
  allInds <- append(allInds, inds)
  betaT_Cfull[inds, ] <- betaT_C
  mu_Cfull[inds, ] <- mu_C
  sigmaE_Cfull[inds, ] <- sigmaE_C
  sigmaP_Cfull[inds, ] <- sigmaP_C
  sigmaPT_Cfull[inds, ] <- sigmaPT_C
  sigmaT_Cfull[inds, ] <- sigmaT_C
}
save(mu_Cfull, betaT_Cfull, sigmaP_Cfull, sigmaT_Cfull,
     sigmaPT_Cfull, sigmaE_Cfull, file = "StanCfullResults.Rdata")

### LOADING FOR GENE SCORING RESULTS

# Set working directory to normal data wd then folder ResultsC
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data/GeneScores")

geneScoresFull <- data.frame(a=numeric(27383), b=numeric(27383), c=numeric(27383))
colnames(geneScoresFull) <- c("logP/T","logP/PT","logPT/T")

nGenes <- 1000
filesPresent <- (1:28)
fileNames <- dir(pattern = "ScoreResults_*")
allInds <- numeric()
for (i in filesPresent) {
  print(paste("Chunk:",i))
  load(fileNames[fileNames == paste("ScoreResults_",i,".Rdata", sep="")])
  inds <- (1:nGenes) + (i - 1) * nGenes
  if (i==28) {
    inds <- inds[1:383]
  }
  allInds <- append(allInds, inds)
  geneScoresFull[inds,] <- geneScores
}
save(geneScoresFull, file = "ScoresFull.Rdata")