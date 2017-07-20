# mark individual genes on histogram of Stan results

library(lme4)
library(arm)
library(rstan)
library(coda)
library(gtools)
library(ggplot2)
library(bayesplot)
library(reshape2)

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
setwd("D:/DataPlus2017/Data")

# load annotation, Stan results, and gene scores
load("myFA.Rdata")
load("StanCfullResults.Rdata")
load("ScoresFull.Rdata")


### INITIALIZE

# # create list of unique gene names and gene regions at each site
# geneNames <- list(866836)
# geneRegions <- list(866836)
# for (i in 1:866836) {
#   if (!(i %% 86684)) {
#     print(paste(i/8668.4,"% complete"))
#   }
#   geneNames[[i]] <- unique(unlist(strsplit(FullAnnotation$UCSC_RefGene_Name[i],split=";")))
#   worklist <- unique(unlist(strsplit(FullAnnotation$UCSC_RefGene_Group[i],split=";")))
#   worklist[worklist=="3'UTR"] <- "3UTR"
#   worklist[worklist=="5'UTR"] <- "5UTR"
#   worklist[worklist=="5URT"] <- "5UTR"
#   if (length(worklist)==0) { worklist <- "blank"}
#   geneRegions[[i]] <- worklist
# }
# save(geneNames, geneRegions, file = "GeneInfo.Rdata")
load("GeneInfo.Rdata")

# make vector of unique genes and regions
uniqueGenes <- unique(unlist(geneNames))
regionTypes <- unique(unlist(geneRegions))

# make dataframes for median and mean values of each parameter
dfMedians <- cbind.data.frame(mu_Cfull$p50, betaT_Cfull$p50, sigmaE_Cfull$p50, sigmaP_Cfull$p50, sigmaPT_Cfull$p50, sigmaT_Cfull$p50)
colnames(dfMedians) <- c("mu","betaT","sigmaE","sigmaP","sigmaPT","sigmaT")
dfMeans <- cbind.data.frame(mu_Cfull$mean, betaT_Cfull$mean, sigmaE_Cfull$mean, sigmaP_Cfull$mean, sigmaPT_Cfull$mean, sigmaT_Cfull$mean)
colnames(dfMeans) <- c("mu","betaT","sigmaE","sigmaP","sigmaPT","sigmaT")

# calculate log ratios: P/T, P/PT, PT/T
logRatios <- cbind.data.frame(log(dfMedians$sigmaP/dfMedians$sigmaT), log(dfMedians$sigmaP/dfMedians$sigmaPT), log(dfMedians$sigmaPT/dfMedians$sigmaT))
colnames(logRatios) <- c("logP/T","logP/PT","logPT/T")
logRatioMeans <- colMeans(logRatios)
# insert means by region

# create histogram plots for each log ratio, to reference later
pPT <- ggplot(logRatios, aes(x=`logP/T`)) + geom_histogram(bins = 100) + ggtitle("Histogram of logP/T") + xlim(c(-3,3)) + ylim(c(0,31000))
pPPT <- ggplot(logRatios, aes(x=`logP/PT`)) + geom_histogram(bins = 100) + ggtitle("Histogram of logP/PT") + xlim(c(-3,3)) + ylim(c(0,31000))
pPTT <- ggplot(logRatios, aes(x=`logPT/T`)) + geom_histogram(bins = 100) + ggtitle("Histogram of logPT/T") + xlim(c(-3,3)) + ylim(c(0,31000))


### FUNCTIONS

inf2NA <- function(x) { x[is.infinite(x)] <- NA; x }

# find indices corresponding to geneString in the list geneNames
geneInd <- function (geneStr) {
  work <- grep(geneStr, geneNames)
  work2 <- sapply(work, function(i) { geneStr %in% geneNames[[i]] } )
  return(work[work2])
}
geneInd <- function (geneStr) {
  return(grep(paste("\\b",geneStr,"\\b",sep=""), geneNames))
}

# plot each site along all three log ratio plots
plotLogRatios <- function (geneStr) {
  gInd <- geneInd(geneStr)
  for (i in colnames(logRatios)){
    if (i=="logP/T") { p1<-pPT }
    if (i=="logP/PT") { p1<-pPPT }
    if (i=="logPT/T") { p1<-pPTT }
    p1 <- p1 + annotate("text", x=logRatios[gInd,i], y = rep(20000,length(gInd)), label="x") + annotate("text", x=0, y=25000, label=(paste(geneStr,"mean:",format(mean(logRatios[gInd,i]),digits=5))))
    print(p1)
  }
}

# returns the functional regions of each index given
plotLogRatiosByRegions <- function (geneStr) {
  gInd <- geneInd(geneStr)
  for (i in colnames(logRatios)){
    if (i=="logP/T") { p1<-pPT }
    if (i=="logP/PT") { p1<-pPPT }
    if (i=="logPT/T") { p1<-pPTT }
    for (j in 1:length(regionTypes)) { # go through each of the region types
      rInd <- integer(0)
      for (k in 1:length(gInd)) { # check for region j at each gInd index k
        if (regionTypes[j] %in% geneRegions[[gInd[k]]]) {
          rInd <- append(rInd,gInd[k])
        }
      }
      col1 <- c("red","orange","gold","green","blue","tomato","sienna","darkviolet")[j]
      p1 <- p1 + annotate("text", x=logRatios[rInd,i], y = rep(j*3000,length(rInd)), label="x", colour=col1) + annotate("text", x=-2, y=j*3000+1000, label=(paste(regionTypes[j],"mean:",format(mean(logRatios[rInd,i]),digits=5))),colour=col1) + annotate("text", x=2.6, y=j*3000+1000, label=(length(rInd)),colour=col1)
    }
    p1 <- p1 + annotate("text", x=logRatios[gInd,i], y = rep(27000,length(gInd)), label=".") + annotate("text", x=0, y=28000, label=(paste(geneStr,"mean:",format(mean(logRatios[gInd,i]),digits=5)))) + annotate("text", x=2, y=29500, label=(paste("Total sites:",length(gInd))))

    print(p1)
  }
}
plotLogRatiosByRegions("GAPDH")

# score gene
scoreGene <- function(geneStr) {
  gInd <- geneInd(geneStr)
  score <- logRatioMeans*0
  for (i in colnames(logRatios)){
    score[i] <- sum(logRatios[gInd,i]-logRatioMeans[i]) / sd(logRatios[gInd,i])
  }
  return(score)
}

# score gene by function region
scoreGeneByRegion <- function(geneStr) {
  
}


### CODE

## Plots

# plot sigmas: E, P, PT, T
df1 <- dfMedians[c("sigmaE","sigmaP","sigmaPT","sigmaT")]
gg <- melt(df1)
ggplot(gg, aes(x=value, fill=variable)) +
  geom_histogram(binwidth=0.01)+xlim(c(0,2))+
  facet_grid(variable~.)

# plot log ratios
gg <- melt(logRatios)
ggplot(gg, aes(x=value, fill=variable)) +
  geom_histogram(binwidth=0.05)+
  facet_grid(variable~.)


## Check individual genes with known conservation

plotLogRatiosByRegions("APC")
plotLogRatiosByRegions("TP53")
plotLogRatiosByRegions("TTN")
plotLogRatiosByRegions("B2M")
plotLogRatiosByRegions("HLA-A")
plotLogRatiosByRegions("HLA-B")

geneCount <- data.frame(table(unlist(geneNames)))
singleGenes <- geneCount$Var1[which(geneCount$Freq == 1)]

# analyze scores with PCA
print(paste("Sites NA:",sum(is.na(geneScoresFull[,1]))))
geneScoresClear <- geneScoresFull[!is.na(geneScoresFull[,1]),]
PCAscore <- prcomp(geneScoresClear)
plot(PCAscore$x[,1],PCAscore$x[,2]) # make a scatterplot
