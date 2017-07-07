# mark individual genes on histogram of Stan results

library(lme4)
library(arm)
library(rstan)
library(coda)
library(gtools)
library(ggplot2)
library(bayesplot)

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
setwd("D:/DataPlus2017/Data")

load("StanCRuns/StanCResults_1.Rdata")

# need stan results first for histogram

load("myFA.Rdata")
# load(stan results at 1000 sites)
load("betaT_C.Rdata")
load("mu_C.Rdata")
load("sigmaE_C.Rdata")
load("sigmaP_C.Rdata")
load("sigmaPT_C.Rdata")
load("sigmaT_C.Rdata")

### Functions for laoding site and Stan fits

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

stanfit3 <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])
  
  
  # Using Model 3: add intra-tumoral variances
  stanFit3 <- stan(file="model3.stan", data=stanDat, control = list(adapt_delta = 0.999))
  
  return(stanFit3=stanFit3)
}

stanfit4 <- function (dataset) {
  
  data <- dataset[dataset[,3]=="T",]
  
  stanDat <- list(pID = as.integer(factor(data$patient)),
                  N = nrow(data),
                  P = nlevels(data$patient),
                  y = data[,1])
  
  
  stanFit4 <- stan(file="model4.stan", data=stanDat, control = list(adapt_delta = 0.999))
  
  return(stanFit4=stanFit4)
}


# find sites with unique names
geneNames <- list(866836)
for (i in 1:866836) {
  geneNames[[i]] <- unique(unlist(strsplit(FullAnnotation$UCSC_RefGene_Name[i],split=";")))
}
# return indices of sites specifically named "APC", and so on
APCind <- grep('APC', geneNames)[(grep('APC', geneNames, value=TRUE)=="APC")]
TP53ind <- grep('TP53', geneNames)[(grep('TP53', geneNames, value=TRUE)=="TP53")]
TTNind <- grep('TTN', geneNames)[(grep('TTN', geneNames, value=TRUE)=="TTN")]
B2Mind <- grep('B2M', geneNames)[(grep('B2M', geneNames, value=TRUE)=="B2M")]
HLAAind <- grep('HLA-A', geneNames)[(grep('HLA-A', geneNames, value=TRUE)=="HLA-A")]
HLABind <- grep('HLA-B', geneNames)[(grep('HLA-B', geneNames, value=TRUE)=="HLA-B")]


## run the indices in stan, store results in logPTratio
# dfC <- data.frame(betaT = numeric(866836), mu = numeric(866836), sigmaE = numeric(866836), sigmaP = numeric(866836), sigmaPT = numeric(866836), sigmaT = numeric(866836), logPTratio = numeric(866836))
# 
# for (i in APCind){
#   data <- site(i)
#   stanC <- stanfit3(data)
#   dfC[i,c("betaT","mu","sigmaE","sigmaP","sigmaPT","sigmaT")] <- summary(stanC)$summary[71:76,1]
#   dfC$logPTratio[i] <- log(dfC$sigmaP[i]/dfC$sigmaT[i])
# }
# 
# for (i in TP53ind){
#   data <- site(i)
#   stanC <- stanfit3(data)
#   dfC[i,c("betaT","mu","sigmaE","sigmaP","sigmaPT","sigmaT")] <- summary(stanC)$summary[71:76,1]
#   dfC$logPTratio[i] <- log(dfC$sigmaP[i]/dfC$sigmaT[i])
# }
# 
# for (i in TTNind){
#   data <- site(i)
#   stanC <- stanfit3(data)
#   dfC[i,c("betaT","mu","sigmaE","sigmaP","sigmaPT","sigmaT")] <- summary(stanC)$summary[71:76,1]
#   dfC$logPTratio[i] <- log(dfC$sigmaP[i]/dfC$sigmaT[i])
# }
# 
# for (i in B2Mind){
#   data <- site(i)
#   stanC <- stanfit3(data)
#   dfC[i,c("betaT","mu","sigmaE","sigmaP","sigmaPT","sigmaT")] <- summary(stanC)$summary[71:76,1]
#   dfC$logPTratio[i] <- log(dfC$sigmaP[i]/dfC$sigmaT[i])
# }
# 
# for (i in HLAAind){
#   data <- site(i)
#   stanC <- stanfit3(data)
#   dfC[i,c("betaT","mu","sigmaE","sigmaP","sigmaPT","sigmaT")] <- summary(stanC)$summary[71:76,1]
#   dfC$logPTratio[i] <- log(dfC$sigmaP[i]/dfC$sigmaT[i])
# }
# 
# for (i in HLABind){
#   data <- site(i)
#   stanC <- stanfit3(data)
#   dfC[i,c("betaT","mu","sigmaE","sigmaP","sigmaPT","sigmaT")] <- summary(stanC)$summary[71:76,1]
#   dfC$logPTratio[i] <- log(dfC$sigmaP[i]/dfC$sigmaT[i])
# }

# save data for later
save(dfS,dfC,file="myDFs.Rdata")
load("myDFs.Rdata")
dfC$logPTTratio <- log(dfC$sigmaPT/dfC$sigmaT)
dfC$logPPTratio <- log(dfC$sigmaP/dfC$sigmaPT)


inf2NA <- function(x) { x[is.infinite(x)] <- NA; x }
# plot histogram of PT ratio, simple and complex model
# first simple model, then complex, for each gene
df <- data.frame(logPTratio = numeric(1000),logPTTratio = numeric(1000),logPPTratio = numeric(1000))
df$logPTratio <- inf2NA(as.matrix(log(sigmaP_C[,1]/sigmaT_C[,1])))
df$logPTTratio <- inf2NA(as.matrix(log(sigmaPT_C[,1]/sigmaT_C[,1])))
df$logPPTratio <- inf2NA(as.matrix(log(sigmaP_C[,1]/sigmaPT_C[,1])))

# APC
ggplot(df, aes(x=logPTratio)) + geom_histogram(bins = 100) + ggtitle("APC locations on P/T ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTratio[APCind], y = rep(c(50,100,150,200,250),11), label="V") + annotate("text", x=0, y= 400, label=(paste("APC mean: ",mean(dfC$logPTratio[APCind], na.rm=TRUE))))
ggplot(df, aes(x=logPTTratio)) + geom_histogram(bins = 100) + ggtitle("APC locations on PT/T ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTTratio[APCind], y = rep(c(50,100,150,200,250),11), label="V") + annotate("text", x=0, y= 400, label=(paste("APC mean: ",mean(dfC$logPTTratio[APCind], na.rm=TRUE))))
ggplot(df, aes(x=logPPTratio)) + geom_histogram(bins = 100) + ggtitle("APC locations on P/PT ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPPTratio[APCind], y = rep(c(50,100,150,200,250),11), label="V") + annotate("text", x=0, y= 400, label=(paste("APC mean: ",mean(dfC$logPPTratio[APCind], na.rm=TRUE))))
# TP53
ggplot(df, aes(x=logPTratio)) + geom_histogram(bins = 100) + ggtitle("TP53 locations on P/T ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTratio[TP53ind], y = rep(c(50,100,150,200,250),3), label="V") + annotate("text", x=0, y= 400, label=(paste("TP53 mean: ",mean(dfC$logPTratio[TP53ind], na.rm=TRUE))))
ggplot(df, aes(x=logPTTratio)) + geom_histogram(bins = 100) + ggtitle("TP53 locations on PT/T ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTTratio[TP53ind], y = rep(c(50,100,150,200,250),3), label="V") + annotate("text", x=0, y= 400, label=(paste("TP53 mean: ",mean(dfC$logPTTratio[TP53ind], na.rm=TRUE))))
ggplot(df, aes(x=logPPTratio)) + geom_histogram(bins = 100) + ggtitle("TP53 locations on P/PT ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPPTratio[TP53ind], y = rep(c(50,100,150,200,250),3), label="V") + annotate("text", x=0, y= 400, label=(paste("TP53 mean: ",mean(dfC$logPPTratio[TP53ind], na.rm=TRUE))))
# TTN
ggplot(df, aes(x=logPTratio)) + geom_histogram(bins = 100) + ggtitle("TTN locations on P/T ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTratio[TTNind], y = rep(c(50,100,150,200,250),7), label="V") + annotate("text", x=0, y= 400, label=(paste("TTN mean: ",mean(dfC$logPTratio[TTNind], na.rm=TRUE))))
ggplot(df, aes(x=logPTTratio)) + geom_histogram(bins = 100) + ggtitle("TTN locations on PT/T ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTTratio[TTNind], y = rep(c(50,100,150,200,250),7), label="V") + annotate("text", x=0, y= 400, label=(paste("TTN mean: ",mean(dfC$logPTTratio[TTNind], na.rm=TRUE))))
ggplot(df, aes(x=logPPTratio)) + geom_histogram(bins = 100) + ggtitle("TTN locations on P/PT ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPPTratio[TTNind], y = rep(c(50,100,150,200,250),7), label="V") + annotate("text", x=0, y= 400, label=(paste("TTN mean: ",mean(dfC$logPPTratio[TTNind], na.rm=TRUE))))
# B2M
ggplot(df, aes(x=logPTratio)) + geom_histogram(bins = 100) + ggtitle("B2M locations on P/T ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTratio[B2Mind], y = rep(c(50,100,150,200,250),5)[-1], label="V") + annotate("text", x=0, y= 400, label=(paste("B2M mean: ",mean(dfC$logPTratio[B2Mind], na.rm=TRUE))))
ggplot(df, aes(x=logPTTratio)) + geom_histogram(bins = 100) + ggtitle("B2M locations on PT/T ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTTratio[B2Mind], y = rep(c(50,100,150,200,250),5)[-1], label="V") + annotate("text", x=0, y= 400, label=(paste("B2M mean: ",mean(dfC$logPTTratio[B2Mind], na.rm=TRUE))))
ggplot(df, aes(x=logPPTratio)) + geom_histogram(bins = 100) + ggtitle("B2M locations on P/PT ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPPTratio[B2Mind], y = rep(c(50,100,150,200,250),5)[-1], label="V") + annotate("text", x=0, y= 400, label=(paste("B2M mean: ",mean(dfC$logPPTratio[B2Mind], na.rm=TRUE))))
# HLA-A
ggplot(df, aes(x=logPTratio)) + geom_histogram(bins = 100) + ggtitle("HLA-A locations on P/T ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTratio[HLAAind], y = rep(c(50,100,150,200,250),7)[-1], label="V") + annotate("text", x=0, y= 400, label=(paste("HLA-A mean: ",mean(dfC$logPTratio[HLAAind], na.rm=TRUE))))
ggplot(df, aes(x=logPTTratio)) + geom_histogram(bins = 100) + ggtitle("HLA-A locations on PT/T ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTTratio[HLAAind], y = rep(c(50,100,150,200,250),7)[-1], label="V") + annotate("text", x=0, y= 400, label=(paste("HLA-A mean: ",mean(dfC$logPTTratio[HLAAind], na.rm=TRUE))))
ggplot(df, aes(x=logPPTratio)) + geom_histogram(bins = 100) + ggtitle("HLA-A locations on P/PT ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPPTratio[HLAAind], y = rep(c(50,100,150,200,250),7)[-1], label="V") + annotate("text", x=0, y= 400, label=(paste("HLA-A mean: ",mean(dfC$logPPTratio[HLAAind], na.rm=TRUE))))
# HLA-B
ggplot(df, aes(x=logPTratio)) + geom_histogram(bins = 100) + ggtitle("HLA-B locations on P/T ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTratio[HLABind], y = rep(c(50,100,150,200,250),9)[-1:-2], label="V") + annotate("text", x=0, y= 400, label=(paste("HLA-B mean: ",mean(dfC$logPTratio[HLABind], na.rm=TRUE))))
ggplot(df, aes(x=logPTTratio)) + geom_histogram(bins = 100) + ggtitle("HLA-B locations on PT/T ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTTratio[HLABind], y = rep(c(50,100,150,200,250),9)[-1:-2], label="V") + annotate("text", x=0, y= 400, label=(paste("HLA-B mean: ",mean(dfC$logPTTratio[HLABind], na.rm=TRUE))))
ggplot(df, aes(x=logPPTratio)) + geom_histogram(bins = 100) + ggtitle("HLA-B locations on P/PT ratio") + xlim(c(-3,3)) + annotate("text", x=dfC$logPPTratio[HLABind], y = rep(c(50,100,150,200,250),9)[-1:-2], label="V") + annotate("text", x=0, y= 400, label=(paste("HLA-B mean: ",mean(dfC$logPPTratio[HLABind], na.rm=TRUE))))







# find sites with unique names
regionNames <- list(866836)
for (i in 1:866836) {
  regionNames[[i]] <- unique(unlist(strsplit(FullAnnotation$UCSC_RefGene_Group[i],split=";")))
}
regionNames[[regionNames=="5URT"]] <- "5'UTR"
# return indices of sites specifically group "TSS200", and so on
TSS200ind <- grep('TSS200', regionNames)[(grep('TSS200', regionNames, value=TRUE)=="TSS200")]
TSS1500ind <- grep('TSS1500', regionNames)[(grep('TSS1500', regionNames, value=TRUE)=="TSS1500")]
UTR3ind <- grep('3UTR', regionNames)[(grep('3UTR', regionNames, value=TRUE)=="3UTR")]
UTR5ind <- grep('5URT', regionNames)[(grep('5URT', regionNames, value=TRUE)=="5URT")]
Bodyind <- grep('Body', regionNames)[(grep('Body', regionNames, value=TRUE)=="Body")]

# plot for APC
df$mu <- mu_C[,1]
APCind2 <- intersect(APCind,Bodyind)
ggplot(df, aes(x=mu)) + geom_histogram(bins = 100) + ggtitle("APC Body locations") + xlim(c(-5,5)) + annotate("text", x=dfC$mu[APCind2], y = rep(30,17), label="V") + annotate("text", x=0, y= 40, label=(paste("APC Body mean: ",mean(dfC$mu[APCind2], na.rm=TRUE))))

