# script for calculating scores on individual genes (27383 total genes)
# using full-scale results from Stan complex model (model3)
# input taskID (1-28) corresponding to 1000 gene chunks starting at site (N-1)*1000+1
# incorporating parallel processing to run through all genes faster
# just need "StanCfullResults.Rdata" in same folder

library(doParallel)

### Initialize

print(paste("Cores: ", detectCores()))

args <- commandArgs(trailingOnly = TRUE)
print(args)
N <- as.numeric(args[1])
print(paste("Task #: ", N))
out.file <- args[2]
log.file <- paste("logTask",N,".txt",sep="")
writeLines(c(""), log.file)

# load data: Stan parameters, geneNames and geneRegions lists
load("StanCfullResults.Rdata")


# make dataframes for median values of each parameter
dfMedians <- cbind.data.frame(mu_Cfull$p50, betaT_Cfull$p50, sigmaE_Cfull$p50, sigmaP_Cfull$p50, sigmaPT_Cfull$p50, sigmaT_Cfull$p50)
colnames(dfMedians) <- c("mu","betaT","sigmaE","sigmaP","sigmaPT","sigmaT")

# calculate log ratios: P/T, P/PT, PT/T
logRatios <- cbind.data.frame(log(dfMedians$sigmaP/dfMedians$sigmaT), log(dfMedians$sigmaP/dfMedians$sigmaPT), log(dfMedians$sigmaPT/dfMedians$sigmaT))
colnames(logRatios) <- c("logP/T","logP/PT","logPT/T")
logRatioMeans <- colMeans(logRatios)
# insert means by region

# Choose 1000 sites by N, and select chunk within datasets
siteInds <- (1:1000) + (N - 1) * 1000
if (N > 27) {
  siteInds <- siteInds[1:383]
}
nsites <- length(siteInds)

print(paste("Making chunk of nsites =", nsites))
chunk <- FullAnnotation[siteInds,]
rm(FullAnnotation)
gc()

### Functions

# find indices corresponding to geneString in the list geneNames
geneInd <- function (geneStr) {
  return(grep(geneStr, geneNames)[(grep(geneStr, geneNames, value=TRUE)==geneStr)])
}

# score gene
scoreGene <- function(geneStr) {
  gInd <- geneInd(geneStr)
  score <- logRatioMeans*0
  for (i in colnames(logRatios)){
    score[i] <- sum(logRatios[gInd,i]-logRatioMeans[i]) / sd(logRatios[gInd,i])
  }
  return(score)
}


### Code


# score all genes
numGenes <- length(uniqueGenes)
geneScores <- data.frame(PTscore=numeric(numGenes), PPTscore=numeric(numGenes), PTTscore=numeric(numGenes))
for (g in 1:numGenes) { # do not run, too slow
  if (!(g %% 1000)) { print(paste("Scoring:",g,"/",numGenes)) }
  geneScores[g,] <- scoreGene(uniqueGenes[g])
}


