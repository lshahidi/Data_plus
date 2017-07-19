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

# load annotation and Stan results
load("myFA.Rdata")
load("StanCfullResults.Rdata")


### INITIALIZE

# create list of unique gene names and gene regions at each site
geneNames <- list(866836)
geneRegions <- list(866836)
for (i in 1:866836) {
  if (!(i %% 86684)) {
    print(paste(i/8668.4,"% complete"))
  }
  geneNames[[i]] <- unique(unlist(strsplit(FullAnnotation$UCSC_RefGene_Name[i],split=";")))
  worklist <- unique(unlist(strsplit(FullAnnotation$UCSC_RefGene_Group[i],split=";")))
  worklist[worklist=="3'UTR"] <- "3UTR"
  worklist[worklist=="5'UTR"] <- "5UTR"
  worklist[worklist=="5URT"] <- "5UTR"
  if (length(worklist)==0) { worklist <- "blank"}
  geneRegions[[i]] <- worklist
}
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
  return(grep(geneStr, geneNames)[(grep(geneStr, geneNames, value=TRUE)==geneStr)])
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
      p1 <- p1 + annotate("text", x=logRatios[rInd,i], y = rep(j*3000,length(rInd)), label="x", colour=col1) + annotate("text", x=-2, y=j*3000+1000, label=(paste(regionTypes[j],"mean:",format(mean(logRatios[rInd,i]),digits=5))),colour=col1)
    }
    p1 <- p1 + annotate("text", x=logRatios[gInd,i], y = rep(27000,length(gInd)), label=".") + annotate("text", x=0, y=28000, label=(paste(geneStr,"mean:",format(mean(logRatios[gInd,i]),digits=5))))

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


### PLOTS

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


# APC
plotLogRatiosByRegions("APC")
scoreGene("APC")
# TP53
plotLogRatiosByRegions("TP53")
scoreGene("TP53")
# TTN
plotLogRatiosByRegions("TTN")
scoreGene("TTN")
# B2M
plotLogRatiosByRegions("B2M")
scoreGene("B2M")
# HLA-A
plotLogRatiosByRegions("HLA-A")
scoreGene("HLA-A")
# HLA-B
plotLogRatiosByRegions("HLA-B")
scoreGene("HLA-B")

# score all genes
numGenes <- length(uniqueGenes)
geneScores <- data.frame(PTscore=numeric(numGenes), PPTscore=numeric(numGenes), PTTscore=numeric(numGenes))
for (g in 1:numGenes) { # do not run, too slow
  if (!(g %% 1000)) { print(paste("Scoring:",g,"/",numGenes)) }
  geneScores[g,] <- scoreGene(uniqueGenes[g])
}
