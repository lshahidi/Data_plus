# script for calculating scores on individual genes (27383 total genes)
# using full-scale results from Stan complex model (model3)
# input taskID (1-28) corresponding to 1000 gene chunks starting at site (N-1)*1000+1
# incorporating parallel processing to run through all genes faster
# just need "StanCfullResults.Rdata" and "GeneInfo" in working directory

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
load("GeneInfo.Rdata")

# make vector of unique genes and regions
uniqueGenes <- unique(unlist(geneNames))
regionTypes <- unique(unlist(geneRegions))

# make dataframes for median values of each parameter
dfMedians <- cbind.data.frame(mu_Cfull$p50, betaT_Cfull$p50, sigmaE_Cfull$p50, sigmaP_Cfull$p50, sigmaPT_Cfull$p50, sigmaT_Cfull$p50)
colnames(dfMedians) <- c("mu","betaT","sigmaE","sigmaP","sigmaPT","sigmaT")

# calculate log ratios: P/T, P/PT, PT/T
logRatios <- cbind.data.frame(log(dfMedians$sigmaP/dfMedians$sigmaT), log(dfMedians$sigmaP/dfMedians$sigmaPT), log(dfMedians$sigmaPT/dfMedians$sigmaT))
colnames(logRatios) <- c("logP/T","logP/PT","logPT/T")
logRatioMeans <- colMeans(logRatios)
# insert means by region

# Choose 1000 genes from uniqueGenes by N, and select geneChunk within gene list
geneIDs <- (1:1000) + (N - 1) * 1000
if (N > 27) {
  geneIDs <- geneIDs[1:383]
}
numGenes <- length(geneIDs)

print(paste("Making chunk of numGenes =", numGenes))
geneChunk <- uniqueGenes[geneIDs]


### Functions

# find indices corresponding to geneString in the list geneNames
geneInd <- function (geneStr) {
  work <- grep(geneStr, geneNames)
  work2 <- sapply(work, function(i) { geneStr %in% geneNames[[i]] } )
  return(work[work2])
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

## score all genes
# run in parallel via doParallel
cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)
print(paste("Cores registered:",getDoParWorkers()))
print(paste("Backend type:",getDoParName()))
print("Starting foreach loop")
ptm <- proc.time()
parData <- foreach(g=iter(1:numGenes), .combine=rbind) %dopar% {
  sink(log.file, append=TRUE)
  print(paste("Gene:", g,"/",numGenes))
  sink()
  scoreGene(geneChunk[g])
}
proc.time() - ptm
stopCluster(cl)

geneScores <- parData
rownames(geneScores) <- geneChunk
rm(parData)
gc()

print(paste("Completed run, now saving"))
# save data
save(geneScores, file = out.file)
