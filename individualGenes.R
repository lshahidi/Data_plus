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


# need stan results first for histogram

load("myFA.Rdata")
# load(stan results at 1000 sites)
load("mu_S.Rdata")
load("sigmaE_S.Rdata")
load("sigmaP_S.Rdata")
load("sigmaE_S.Rdata")
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


# run the indices in stan, store results in logPTratio
dfS <- data.frame(mu = numeric(866836), sigmaE = numeric(866836), sigmaP = numeric(866836), sigmaT = numeric(866836), logPTratio = numeric(866836))
dfC <- data.frame(betaT = numeric(866836), mu = numeric(866836), sigmaE = numeric(866836), sigmaP = numeric(866836), sigmaPT = numeric(866836), sigmaT = numeric(866836), logPTratio = numeric(866836))

for (i in APCind){
  data <- site(i)
  stanS <- stanfit4(data)
  dfS[i,c("mu","sigmaE","sigmaP","sigmaT")] <- summary(stanS)$summary[49:52,1]
  dfS$logPTratio[i] <- log(dfS$sigmaP[i]/dfS$sigmaT[i])
  stanC <- stanfit3(data)
  dfC[i,c("betaT","mu","sigmaE","sigmaP","sigmaPT","sigmaT")] <- summary(stanC)$summary[71:76,1]
  dfC$logPTratio[i] <- log(dfC$sigmaP[i]/dfC$sigmaT[i])
}

for (i in TP53ind){
  data <- site(i)
  stanS <- stanfit4(data)
  dfS[i,c("mu","sigmaE","sigmaP","sigmaT")] <- summary(stanS)$summary[49:52,1]
  dfS$logPTratio[i] <- log(dfS$sigmaP[i]/dfS$sigmaT[i])
  stanC <- stanfit3(data)
  dfC[i,c("betaT","mu","sigmaE","sigmaP","sigmaPT","sigmaT")] <- summary(stanC)$summary[71:76,1]
  dfC$logPTratio[i] <- log(dfC$sigmaP[i]/dfC$sigmaT[i])
}

for (i in TTNind){
  data <- site(i)
  stanS <- stanfit4(data)
  dfS[i,c("mu","sigmaE","sigmaP","sigmaT")] <- summary(stanS)$summary[49:52,1]
  dfS$logPTratio[i] <- log(dfS$sigmaP[i]/dfS$sigmaT[i])
  stanC <- stanfit3(data)
  dfC[i,c("betaT","mu","sigmaE","sigmaP","sigmaPT","sigmaT")] <- summary(stanC)$summary[71:76,1]
  dfC$logPTratio[i] <- log(dfC$sigmaP[i]/dfC$sigmaT[i])
}

for (i in B2Mind){
  data <- site(i)
  stanS <- stanfit4(data)
  dfS[i,c("mu","sigmaE","sigmaP","sigmaT")] <- summary(stanS)$summary[49:52,1]
  dfS$logPTratio[i] <- log(dfS$sigmaP[i]/dfS$sigmaT[i])
  stanC <- stanfit3(data)
  dfC[i,c("betaT","mu","sigmaE","sigmaP","sigmaPT","sigmaT")] <- summary(stanC)$summary[71:76,1]
  dfC$logPTratio[i] <- log(dfC$sigmaP[i]/dfC$sigmaT[i])
}

for (i in HLAAind){
  data <- site(i)
  stanS <- stanfit4(data)
  dfS[i,c("mu","sigmaE","sigmaP","sigmaT")] <- summary(stanS)$summary[49:52,1]
  dfS$logPTratio[i] <- log(dfS$sigmaP[i]/dfS$sigmaT[i])
  stanC <- stanfit3(data)
  dfC[i,c("betaT","mu","sigmaE","sigmaP","sigmaPT","sigmaT")] <- summary(stanC)$summary[71:76,1]
  dfC$logPTratio[i] <- log(dfC$sigmaP[i]/dfC$sigmaT[i])
}

for (i in HLABind){
  data <- site(i)
  stanS <- stanfit4(data)
  dfS[i,c("mu","sigmaE","sigmaP","sigmaT")] <- summary(stanS)$summary[49:52,1]
  dfS$logPTratio[i] <- log(dfS$sigmaP[i]/dfS$sigmaT[i])
  stanC <- stanfit3(data)
  dfC[i,c("betaT","mu","sigmaE","sigmaP","sigmaPT","sigmaT")] <- summary(stanC)$summary[71:76,1]
  dfC$logPTratio[i] <- log(dfC$sigmaP[i]/dfC$sigmaT[i])
}

# save data for later
save(dfS,dfC,file="myDFs.Rdata")
load("myDFs.Rdata")


inf2NA <- function(x) { x[is.infinite(x)] <- NA; x }
# plot histogram of PT ratio, simple and complex model
# first simple model, then complex, for each gene
df <- data.frame(logPTratioS1K = numeric(1000), logPTratioC1K = numeric(1000))
df$logPTratioS1K <- inf2NA(as.matrix(log(sigmaP_S[,1]/sigmaT_S[,1])))
df$logPTratioC1K <- inf2NA(as.matrix(log(sigmaP_C[,1]/sigmaT_C[,1])))
# APC
ggplot(df, aes(x=logPTratioS1K)) + geom_histogram(bins = 100) + ggtitle("APC locations") + xlim(c(-3,3)) + annotate("text", x=dfS$logPTratio[APCind], y = rep(c(200,250,300,350,400),11), label="V") + annotate("text", x=0, y= 600, label=(paste("APC mean: ",mean(dfS$logPTratio[APCind], na.rm=TRUE))))
ggplot(df, aes(x=logPTratioC1K)) + geom_histogram(bins = 100) + ggtitle("APC locations") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTratio[APCind], y = rep(c(50,100,150,200,250),11), label="V") + annotate("text", x=0, y= 400, label=(paste("APC mean: ",mean(dfC$logPTratio[APCind], na.rm=TRUE))))
# TP53
ggplot(df, aes(x=logPTratioS1K)) + geom_histogram(bins = 100) + ggtitle("TP53 locations") + xlim(c(-3,3)) + annotate("text", x=dfS$logPTratio[TP53ind], y = rep(c(200,250,300,350,400),3), label="V") + annotate("text", x=0, y= 600, label=(paste("TP53 mean: ",mean(dfS$logPTratio[TP53ind], na.rm=TRUE))))
ggplot(df, aes(x=logPTratioC1K)) + geom_histogram(bins = 100) + ggtitle("TP53 locations") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTratio[TP53ind], y = rep(c(50,100,150,200,250),3), label="V") + annotate("text", x=0, y= 400, label=(paste("TP53 mean: ",mean(dfC$logPTratio[TP53ind], na.rm=TRUE))))
# TTN
ggplot(df, aes(x=logPTratioS1K)) + geom_histogram(bins = 100) + ggtitle("TTN locations") + xlim(c(-3,3)) + annotate("text", x=dfS$logPTratio[TTNind], y = rep(c(200,250,300,350,400),7), label="V") + annotate("text", x=0, y= 600, label=(paste("TTN mean: ",mean(dfS$logPTratio[TTNind], na.rm=TRUE))))
ggplot(df, aes(x=logPTratioC1K)) + geom_histogram(bins = 100) + ggtitle("TTN locations") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTratio[TTNind], y = rep(c(50,100,150,200,250),7), label="V") + annotate("text", x=0, y= 400, label=(paste("TTN mean: ",mean(dfC$logPTratio[TTNind], na.rm=TRUE))))
# B2M
ggplot(df, aes(x=logPTratioS1K)) + geom_histogram(bins = 100) + ggtitle("B2M locations") + xlim(c(-3,3)) + annotate("text", x=dfS$logPTratio[B2Mind], y = rep(c(200,250,300,350,400),5)[-1], label="V") + annotate("text", x=0, y= 600, label=(paste("B2M mean: ",mean(dfS$logPTratio[B2Mind], na.rm=TRUE))))
ggplot(df, aes(x=logPTratioC1K)) + geom_histogram(bins = 100) + ggtitle("B2M locations") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTratio[B2Mind], y = rep(c(50,100,150,200,250),5)[-1], label="V") + annotate("text", x=0, y= 400, label=(paste("B2M mean: ",mean(dfC$logPTratio[B2Mind], na.rm=TRUE))))
# HLA-A
ggplot(df, aes(x=logPTratioS1K)) + geom_histogram(bins = 100) + ggtitle("HLA-A locations") + xlim(c(-3,3)) + annotate("text", x=dfS$logPTratio[HLAAind], y = rep(c(200,250,300,350,400),7)[-1], label="V") + annotate("text", x=0, y= 600, label=(paste("HLA-A mean: ",mean(dfS$logPTratio[HLAAind], na.rm=TRUE))))
ggplot(df, aes(x=logPTratioC1K)) + geom_histogram(bins = 100) + ggtitle("HLA-A locations") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTratio[HLAAind], y = rep(c(50,100,150,200,250),7)[-1], label="V") + annotate("text", x=0, y= 400, label=(paste("HLA-A mean: ",mean(dfC$logPTratio[HLAAind], na.rm=TRUE))))
# HLA-B
ggplot(df, aes(x=logPTratioS1K)) + geom_histogram(bins = 100) + ggtitle("HLA-B locations") + xlim(c(-3,3)) + annotate("text", x=dfS$logPTratio[HLABind], y = rep(c(200,250,300,350,400),9)[-1:-2], label="V") + annotate("text", x=0, y= 600, label=(paste("HLA-B mean: ",mean(dfS$logPTratio[HLABind], na.rm=TRUE))))
ggplot(df, aes(x=logPTratioC1K)) + geom_histogram(bins = 100) + ggtitle("HLA-B locations") + xlim(c(-3,3)) + annotate("text", x=dfC$logPTratio[HLABind], y = rep(c(50,100,150,200,250),9)[-1:-2], label="V") + annotate("text", x=0, y= 400, label=(paste("HLA-B mean: ",mean(dfC$logPTratio[HLABind], na.rm=TRUE))))

