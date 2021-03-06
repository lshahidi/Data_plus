# this is the first code for modelling our data with lme4
# will incorporate RStan model as well


### INITIALIZE

library(lme4)
library(arm)
library(rstan)
library(coda)
library(gtools)
library(reshape2)
library(ggplot2)
library(bayesplot)

# here set the working directory that points to the data folder
# e.g. the folder with annotated data saved as "myFA.Rdata"
# please just comment out the directories already here, so each
# user can uncomment the setwd for them as we move code around

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
setwd("D:/DataPlus2017/Data")


### LOAD DATA

# load fully annotated data (saved from LoadDataAndQC.R)
load("myFA.Rdata")
load("myLMERVars.Rdata")

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

fit1 <- lmer(X1 ~ tInd + (tInd|patient), Data1)

# fitGland <- lmer(X1 ~ tInd + (tInd|patient) + (0+gInd|patient:side), Data1)


## extract variance at 1000 sites
# nSites<-1000
# sigmaLMER <- data.frame(sigmaT=numeric(nSites),sigmaP=numeric(nSites),sigmaE=numeric(nSites))
# 
# ptm <- proc.time()
# for (i in 1:nSites){
#   print(paste("site ",i))
#   
#   DataI <- site(i)
#   names(DataI)[1] <- "beta"
#   
#   fit <- lmer(beta ~ tInd + (tInd|patient), DataI)
#   
#   sigmaLMER$sigmaT[i] <- as.data.frame(VarCorr(fit))$sdcor[1]
#   sigmaLMER$sigmaP[i] <- as.data.frame(VarCorr(fit))$sdcor[2]
#   sigmaLMER$sigmaE[i] <- as.data.frame(VarCorr(fit))$sdcor[4]
# }
# proc.time() - ptm


# plot variances
sigmaLMER2 <- sigmaLMER[,c("sigmaP","sigmaT","sigmaE")]
sigmaLMER2[sigmaLMER2>2] <- NA
sigmaLMER2[sigmaLMER2==0] <- NA
gg <- melt(sigmaLMER2)
ggplot(gg, aes(x=value, fill=variable)) +
  geom_histogram(binwidth=0.05)+
  facet_grid(variable~.)

# log sigmaT
logSigT <- log(sigmaLMER2$sigmaT)
hist(logSigT[logSigT>(-30)],100)
# log sigmaP
hist(log(sigmaLMER2$sigmaP),100)

# examine P/T ratio
PTratio <- sigmaLMER$sigmaP/sigmaLMER$sigmaT
PTratio[PTratio>100] <- NA
hist(PTratio,100)

logPTratio <- log(sigmaLMER$sigmaP/sigmaLMER$sigmaT)
logPTratio[logPTratio<(-50)] <- NA
logPTratio[logPTratio>50] <- NA
hist(logPTratio,100)

# layout boxplot is at the bottom 
nf <- layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(3,1))
par(mar=c(3.1, 3.1, 1.1, 2.1))
hist(logPTratio,xlim=c(-15,15),100)
boxplot(logPTratio, horizontal=TRUE,  outline=TRUE, ylim=c(-15,15), frame=F, width = 10)

mean(PTratio, na.rm=TRUE)
mean(log(PTratio), na.rm=TRUE)

sigmaLMER$logPTratio <- logPTratio

inf2NA <- function(x) { x[is.infinite(x)] <- NA; x }
# order by sigmaP
ordSigP <- sigmaLMER[order(sigmaLMER$sigmaP),]
ordSigP$logPTratio <- inf2NA(ordSigP$logPTratio)
# order by sigmaT
ordSigT <- sigmaLMER[order(sigmaLMER$sigmaT),]
ordSigT$logPTratio <- inf2NA(ordSigT$logPTratio)

# annotate conserved sites, with low sigmaP
nTail = 8668 #866836/100
lowP <- head(ordSigP, nTail)
# number of methylated/unmethylated, based on mu [false true]
as.numeric(summary(lowP$mu>=0)[2:3])/8668
as.numeric(summary(ordSigP$mu>=0)[2:3])/866836

# plot variances again for lowP
lowPsigs <- lowP[,c("sigmaP","sigmaT","sigmaE")]
lowPsigs[lowPsigs>2] <- NA
gg <- melt(lowPsigs)
ggplot(gg, aes(x=value, fill=variable)) +
  geom_histogram(binwidth=0.05)+
  facet_grid(variable~.)

# print annotated gene names for lowP sites
lowPnames <- FullAnnotation$UCSC_RefGene_Name[as.numeric(row.names(lowP))]
uniqueNames <- list()
for (i in 1:8668) {
  uniqueNames <- c(uniqueNames, unique(unlist(strsplit(lowPnames[i],split=";"))))
}
uniqueNames <- unique(unlist(uniqueNames))

# order by sigma T, select top 10%
lowPordSigT <- lowP[order(lowP$sigmaT),]
lowPhiT <- tail(lowPordSigT,2000)
lowPhiTnames <- FullAnnotation$UCSC_RefGene_Name[as.numeric(row.names(lowPhiT))]
uniqueNames <- list()
for (i in 1:8668) {
  uniqueNames <- c(uniqueNames, unique(unlist(strsplit(lowPhiTnames[i],split=";"))))
}
uniqueNames <- unique(unlist(uniqueNames))


# sigmaT vs sigmaP vs sigmaE
ggplot(data = ordSigT, aes(x=sigmaT, y=sigmaP)) + geom_point(alpha = 0.05)
ggplot(data = ordSigT, aes(x=sigmaT, y=sigmaE)) + geom_point(alpha = 0.05)
ggplot(data = ordSigT, aes(x=sigmaP, y=sigmaE)) + geom_point(alpha = 0.05)
cor(ordSigT$sigmaT,ordSigT$sigmaP, use = "pairwise.complete.obs")
cor(ordSigT$sigmaT,ordSigT$sigmaE, use = "pairwise.complete.obs")
cor(ordSigT$sigmaP,ordSigT$sigmaE, use = "pairwise.complete.obs")
# vs logPTratio might be useful
ggplot(data = ordSigT, aes(x=sigmaT, y=logPTratio)) + geom_point(alpha = 0.05)

# figure 3 - sort sigma T (done above) and plot in order with other sigmas
ggplot(data = ordSigT, aes(x=seq_along(ordSigT$sigmaT), y=sigmaT)) + geom_point(alpha = 0.05) + labs(x="Index, ordered by sigmaT")
ggplot(data = ordSigT, aes(x=seq_along(ordSigT$sigmaT), y=sigmaP)) + geom_point(alpha = 0.05) + labs(x="Index, ordered by sigmaT")
ggplot(data = ordSigT, aes(x=seq_along(ordSigT$sigmaT), y=sigmaE)) + geom_point(alpha = 0.05) + labs(x="Index, ordered by sigmaT")
ggplot(data = ordSigT, aes(x=seq_along(ordSigT$sigmaT), y=logPTratio)) + geom_point(alpha = 0.05) + labs(x="Index, ordered by sigmaT")

# same for sigma p
ggplot(data = ordSigP, aes(x=seq_along(ordSigP$sigmaP), y=sigmaP)) + geom_point(alpha = 0.05) + labs(x="Index, ordered by sigmaP")
ggplot(data = ordSigP, aes(x=seq_along(ordSigP$sigmaP), y=sigmaT)) + geom_point(alpha = 0.05) + labs(x="Index, ordered by sigmaP")
ggplot(data = ordSigP, aes(x=seq_along(ordSigP$sigmaP), y=sigmaE)) + geom_point(alpha = 0.05) + labs(x="Index, ordered by sigmaP")
ggplot(data = ordSigP, aes(x=seq_along(ordSigP$sigmaP), y=logPTratio)) + geom_point(alpha = 0.05) + labs(x="Index, ordered by sigmaP")


cor(seq_along(ordSigT$sigmaT),ordSigT$sigmaP, use = "pairwise.complete.obs")
cor(seq_along(ordSigT$sigmaT),ordSigT$sigmaE, use = "pairwise.complete.obs")
cor(seq_along(ordSigT$sigmaT),ordSigT$logPTratio, use = "pairwise.complete.obs")
cor(seq_along(ordSigP$sigmaP),ordSigP$sigmaT, use = "pairwise.complete.obs")
cor(seq_along(ordSigP$sigmaP),ordSigP$sigmaE, use = "pairwise.complete.obs")
cor(seq_along(ordSigP$sigmaP),ordSigP$logPTratio, use = "pairwise.complete.obs")


# take all zero sigmaT
zeroSigP <- head(ordSigP,112790)
gg <- melt(zeroSigP[c("sigmaT","sigmaE")])
ggplot(gg, aes(x=value, fill=variable)) +
  geom_histogram(binwidth=0.05)+
  facet_grid(variable~.)
zeroSigPraw <- FullAnnotation[as.numeric(row.names(zeroSigP)),]
zeroSigPraw2 <- site(as.numeric(row.names(zeroSigP)))


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

# plot histogram of PT ratio
df <- data.frame(inf2NA(logPTratio))
ggplot(df, aes(x=logPTratio)) + geom_histogram(bins = 100) + ggtitle("APC locations") + xlim(c(-30,30)) + annotate("text", x=logPTratio[APCind], y = rep(c(20000,25000,30000,35000,40000),11), label="V") + annotate("text", x=8, y= 60000, label=(paste("APC mean: ",mean(logPTratio[APCind], na.rm=TRUE))))
ggplot(df, aes(x=logPTratio)) + geom_histogram(bins = 100) + ggtitle("TP53 locations") + xlim(c(-30,30)) + annotate("text", x=logPTratio[TP53ind], y = rep(c(20000,25000,30000,35000,40000),3), label="V") + annotate("text", x=8, y= 60000, label=(paste("TP53 mean: ",mean(logPTratio[TP53ind], na.rm=TRUE))))
ggplot(df, aes(x=logPTratio)) + geom_histogram(bins = 100) + ggtitle("TTN locations") + xlim(c(-30,30)) + annotate("text", x=logPTratio[TTNind], y = rep(c(20000,25000,30000,35000,40000),7), label="V") + annotate("text", x=8, y= 60000, label=(paste("TTN mean: ",mean(logPTratio[TTNind], na.rm=TRUE))))
ggplot(df, aes(x=logPTratio)) + geom_histogram(bins = 100) + ggtitle("B2M locations") + xlim(c(-30,30)) + annotate("text", x=logPTratio[B2Mind], y = rep(c(20000,25000,30000,35000,40000),5)[-1], label="V") + annotate("text", x=8, y= 60000, label=(paste("B2M mean: ",mean(logPTratio[B2Mind], na.rm=TRUE))))
ggplot(df, aes(x=logPTratio)) + geom_histogram(bins = 100) + ggtitle("HLA-A locations") + xlim(c(-30,30)) + annotate("text", x=logPTratio[HLAAind], y = rep(c(20000,25000,30000,35000,40000),7)[-1], label="V") + annotate("text", x=8, y= 60000, label=(paste("HLA-A mean: ",mean(logPTratio[HLAAind], na.rm=TRUE))))
ggplot(df, aes(x=logPTratio)) + geom_histogram(bins = 100) + ggtitle("HLA-B locations") + xlim(c(-30,30)) + annotate("text", x=logPTratio[HLABind], y = rep(c(20000,25000,30000,35000,40000),9)[-1:-2], label="V") + annotate("text", x=8, y= 60000, label=(paste("HLA-B mean: ",mean(logPTratio[HLABind], na.rm=TRUE))))




# Variances analysis by functional regions (figure 2)

load("myLMERVars.Rdata")
sigmaLMER$logr <- log(sigmaLMER$sigmaP/sigmaLMER$sigmaT)

# identify enhancers
Enhancer <- rep(NA,dim(FullAnnotation)[1])
for (i in 1:dim(FullAnnotation)[1]) {
  if (FullAnnotation$Phantom4_Enhancers[i] != '' | FullAnnotation$Phantom5_Enhancers[i] != '') {
    Enhancer[i] <- 1
  } else {
    Enhancer[i] <- 0
  }
}

# identify all possibilities in gene groups
temp <- strsplit(FullAnnotation$UCSC_RefGene_Group,";")
unique(unlist(temp))

# [1] "TSS1500" "Body"    "3'UTR"   "1stExon" "TSS200"  "5'UTR"   "5URT"    "3UTR"   
# [9] "ExonBnd"

Promoter <- rep(NA,dim(FullAnnotation)[1])
Body <- rep(NA,dim(FullAnnotation)[1])
Exon <- rep(NA,dim(FullAnnotation)[1])
UTR_5 <- rep(NA,dim(FullAnnotation)[1])
UTR_3 <- rep(NA,dim(FullAnnotation)[1])

for (i in 1:dim(FullAnnotation)[1]) {
  split <- strsplit(FullAnnotation$UCSC_RefGene_Group[i],";")
  
  if (sum(split[[1]] %in% "TSS1500" | split[[1]] %in% "TSS200") != 0) {
    Promoter[i] <- 1
  } else {
    Promoter[i] <- 0
  }
  if (sum(split[[1]] %in% "Body") != 0 ) {
    Body[i] <- 1
  } else {
    Body[i] <- 0
  }
  if (sum(split[[1]] %in% "1stExon" | split[[1]] %in% "ExonBnd") != 0 ) {
    Exon[i] <- 1
  } else {
    Exon[i] <- 0
  }
  if (sum(split[[1]] %in% "5URT" | split[[1]] %in% "5'UTR") != 0 ) {
    UTR_5[i] <- 1
  } else {
    UTR_5[i] <- 0
  }
  if (sum(split[[1]] %in% "3UTR" | split[[1]] %in% "3'UTR") != 0 ) {
    UTR_3[i] <- 1
  } else {
    UTR_3[i] <- 0
  }
}

sigmaLMER$Enhancer <- Enhancer
sigmaLMER$Promoter <- Promoter
sigmaLMER$Exon <- Exon
sigmaLMER$Body <- Body
sigmaLMER$UTR_3 <- UTR_3
sigmaLMER$UTR_5 <- UTR_5

# Exclude inf values
sigmaLMER <- sigmaLMER[sigmaLMER$logr != Inf,]

# Seperate plots

par(mfrow=c(2,3))

hist(sigmaLMER$logr[sigmaLMER$Enhancer==1],
     main="Distribution of Enhancers",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=1000,xlim=c(-10,10),ylim=c(0,20000),
     col="cyan")

hist(sigmaLMER$logr[sigmaLMER$Promoter==1],
     main="Distribution of Promoters",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=1000,xlim=c(-10,10),ylim=c(0,20000),
     col="coral")

hist(sigmaLMER$logr[sigmaLMER$Exon==1],
     main="Distribution of Exon Regions",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=1000,xlim=c(-10,10),ylim=c(0,20000),
     col="darkorchid1")

hist(sigmaLMER$logr[sigmaLMER$Body==1],
     main="Distribution of Gene Body",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=1000,xlim=c(-10,10),ylim=c(0,20000),
     col="deeppink")

hist(sigmaLMER$logr[sigmaLMER$UTR_3==1],
     main="Distribution of 3'UTR",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=1000,xlim=c(-10,10),ylim=c(0,15000),
     col="dodgerblue")

hist(sigmaLMER$logr[sigmaLMER$UTR_5==1],
     main="Distribution of 5'UTR",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=1000,xlim=c(-10,10),ylim=c(0,15000),
     col="gold")

# All plots

par(mfrow=c(1,1))

hist(sigmaLMER$logr,main="Distribution of Enhancers",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=1000,xlim=c(-15,15),ylim=c(0,50000),
     col=grey.colors(1,alpha=0.1))

hist(sigmaLMER$logr[sigmaLMER$Enhancer==1],
     breaks=1000,col="cyan", add=TRUE)

legend("topright",c("All Sites","Enhancer"),fill=c(grey.colors(1,alpha=0.1),"cyan"))


hist(sigmaLMER$logr,main="Distribution of Promoters",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=1000,xlim=c(-15,15),ylim=c(0,50000),
     col=grey.colors(1,alpha=0.1))

hist(sigmaLMER$logr[sigmaLMER$Promoter==1],
     breaks=1000,col="coral", add=TRUE)

legend("topright",c("All Sites","Promoter"),fill=c(grey.colors(1,alpha=0.1),"coral"))


hist(sigmaLMER$logr,main="Distribution of Exon Regions",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=1000,xlim=c(-15,15),ylim=c(0,50000),
     col=grey.colors(1,alpha=0.1))

hist(sigmaLMER$logr[sigmaLMER$Exon==1],
     breaks=1000,col="darkorchid1", add=TRUE)

legend("topright",c("All Sites","Exon Region"),fill=c(grey.colors(1,alpha=0.1),"darkorchid1"))


hist(sigmaLMER$logr,main="Distribution of Gene Body",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=1000,xlim=c(-15,15),ylim=c(0,50000),
     col=grey.colors(1,alpha=0.1))

hist(sigmaLMER$logr[sigmaLMER$Body==1],
     breaks=1000,col="deeppink", add=TRUE)

legend("topright",c("All Sites","Gene Body"),fill=c(grey.colors(1,alpha=0.1),"deeppink"))



hist(sigmaLMER$logr,main="Distribution of 3'UTR",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=1000,xlim=c(-15,15),ylim=c(0,50000),
     col=grey.colors(1,alpha=0.1))

hist(sigmaLMER$logr[sigmaLMER$UTR_3==1],
     breaks=1000,col="dodgerblue", add=TRUE)

legend("topright",c("All Sites","3'UTR"),fill=c(grey.colors(1,alpha=0.1),"dodgerblue"))


hist(sigmaLMER$logr,main="Distribution of 5'UTR",
     xlab="log ratio of sigmaP/sigmaT",
     breaks=1000,xlim=c(-15,15),ylim=c(0,50000),
     col=grey.colors(1,alpha=0.1))

hist(sigmaLMER$logr[sigmaLMER$UTR_5==1],
     breaks=1000,col="gold", add=TRUE)

legend("topright",c("All Sites","5'UTR"),fill=c(grey.colors(1,alpha=0.1),"gold"))






# Analysis 1

r_enhancer <- sigmaLMER$logr[sigmaLMER$Enhancer==1]
r_promoter <- sigmaLMER$logr[sigmaLMER$Promoter==1]
r_exon <- sigmaLMER$logr[sigmaLMER$Exon==1]
r_body <- sigmaLMER$logr[sigmaLMER$Body==1]
r_3 <- sigmaLMER$logr[sigmaLMER$UTR_3==1]
r_5 <- sigmaLMER$logr[sigmaLMER$UTR_5==1]


test.data <- data.frame(y=c(r_enhancer,r_promoter,r_exon,r_body,r_3,r_5),
                        region=factor(rep(c("Enhancer","Promoter","Exon","Body",
                                     "3'UTR","5'UTR"),times=c(length(r_enhancer),
                                                              length(r_promoter),
                                                              length(r_exon),
                                                              length(r_body),
                                                              length(r_3),
                                                              length(r_5)))))

# Elimintate infinite values
test.data <- test.data[test.data$y != Inf,]



mean(test.data$y[test.data$region=="Enhancer"])
mean(test.data$y[test.data$region=="Promoter"])
mean(test.data$y[test.data$region=="Exon"])
mean(test.data$y[test.data$region=="Body"])
mean(test.data$y[test.data$region=="3'UTR"])
mean(test.data$y[test.data$region=="5'UTR"])


median(test.data$y[test.data$region=="Enhancer"])
median(test.data$y[test.data$region=="Promoter"])
median(test.data$y[test.data$region=="Exon"])
median(test.data$y[test.data$region=="Body"])
median(test.data$y[test.data$region=="3'UTR"])
median(test.data$y[test.data$region=="5'UTR"])



test <- aov(y ~ region, data=test.data)
anova(test)

# Tukey
TukeyHSD(test)

# Model Checking
par(mfrow=c(2,2))
plot(test)


# Kruskal-Wallis Test
kruskal.test(y~region, data=test.data)

# Dunn test for multiple comparison 
dunn.test(test.data$y,test.data$region, method="bh")






### FIT WITH STAN


# Function to build stan model for each site
stanfit2 <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])
  
  
  # Model 1: uniform priors
  # stanFit <- stan(file="model.stan",data=stanDat)
  
  
  # Using Model 2: wide range uniform prior on sd and flat normal on fixed effects
  stanFit2 <- stan(file="model2.stan", data=stanDat, control = list(adapt_delta = 0.999))
  
  
  return(stanFit2=stanFit2)
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

Data1 <- site(1)
Data2 <- site(2)
Data3 <- site(3)
Data4 <- site(4)
Data5 <- site(5)
Data6 <- site(6)
Data7 <- site(7)
Data8 <- site(8)


stan1 <- stanfit2(Data1)
stan1t <- stanfit3(Data1)

plot(stan1, pars=c("sigma_t","sigma_p","sigma_e"))
plot(stan1t, pars=c("sigma_t","sigma_p","sigma_e","sigma_pt"))





# Explore posterior distributions

posterior1 <- as.array(stan1)
  
mcmc_areas(
    posterior1, 
    pars = c("sigma_e", "sigma_p", "sigma_t"),
    prob = 0.8, # 80% intervals
    prob_outer = 0.95, 
    point_est = "mean"
)
  

posterior2 <- as.array(stan2)

mcmc_areas(
  posterior2, 
  pars = c("sigma_e", "sigma_p", "sigma_t"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, 
  point_est = "mean"
)


posterior3 <- as.array(stan3)

mcmc_areas(
  posterior3, 
  pars = c("sigma_e", "sigma_p", "sigma_t"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, 
  point_est = "mean"
)


posterior4 <- as.array(stan4)

mcmc_areas(
  posterior4, 
  pars = c("sigma_e", "sigma_p", "sigma_t"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, 
  point_est = "mean"
)
  

posterior5 <- as.array(stan5)

mcmc_areas(
  posterior5, 
  pars = c("sigma_e", "sigma_p", "sigma_t"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, 
  point_est = "mean"
)
  


# Extract estimates from stan models

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

summary(stan1)$summary[33:37,1]

getmode(posterior1[1:1000,1:4,33])
getmode(posterior1[1:1000,1:4,37])


summary(stan2)$summary[33:37,1]

rbind(summary(stan1)$summary[33:37,1],summary(stan2)$summary[33:37,1])



# ALL SITES
# Test on first ten sites

nsites <- 10

betaT <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                    p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                    p97.5=numeric(nsites))

mu <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                 p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                 p97.5=numeric(nsites))

sigma_e <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                      p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                      p97.5=numeric(nsites))


sigma_p <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                      p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                      p97.5=numeric(nsites))


sigma_t <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                      p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                      p97.5=numeric(nsites))


for (i in 1:10) {
  data <- site(i)
  stan <- stanfit(data)
  posterior <- as.array(stan)
  
  betaT[i,1] <- summary(stan)$summary[33,1]
  betaT[i,2] <- getmode(posterior[1:1000,1:4,33])
  betaT[i,3:7] <- summary(stan)$summary[33,3:7]
  
  mu[i,1] <- summary(stan)$summary[34,1]
  mu[i,2] <- getmode(posterior[1:1000,1:4,34])
  mu[i,3:7] <- summary(stan)$summary[34,3:7]
  
  sigma_e[i,1] <- summary(stan)$summary[35,1]
  sigma_e[i,2] <- getmode(posterior[1:1000,1:4,35])
  sigma_e[i,3:7] <- summary(stan)$summary[35,3:7]
  
  
  sigma_p[i,1] <- summary(stan)$summary[36,1]
  sigma_p[i,2] <- getmode(posterior1[1:1000,1:4,36])
  sigma_p[i,3:7] <- summary(stan)$summary[36,3:7]
  
  
  sigma_t[i,1] <- summary(stan)$summary[37,1]
  sigma_t[i,2] <- getmode(posterior1[1:1000,1:4,37])
  sigma_t[i,3:7] <- summary(stan)$summary[37,3:7]
}




# Comparing lmer and stan results, based on random 200 sites

# Randomly choose 200 sites

set.seed(123)
indice <- sample(1:dim(FullAnnotation)[1],200,replace=FALSE)
indice <- sort(indice)


# lmer

sigmaLMER_200 <- sigmaLMER[indice,1:3]

# stan

sigmaSTAN_200 <- as.data.frame(matrix(NA,200,3))
count <- 1

for (i in indice) {
  data <- site(i)
  fit <- stanfit(data)
  sigmaSTAN_200[count,] <- summary(fit)$summary[35:37,1]
  
  count <- count+1
}



colnames(sigmaSTAN_200) <- c("sigmaE","sigmaP","sigmaT")
rownames(sigmaSTAN_200) <- indice

sigmaSTAN_200 <- sigmaSTAN_200[,c(2,3,1)]


# Variance sum
sigmaLMER_200$var_sum <- apply(sigmaLMER_200^2,1,sum)
sigmaSTAN_200$var_sum <- apply(sigmaSTAN_200^2,1,sum)

# lmer vs. stan
sigmaDIFF <- sigmaLMER_200 - sigmaSTAN_200


# Difference mean
apply(sigmaDIFF,2,mean)

# Difference sd
sqrt(apply(sigmaDIFF,2,var))


# Plots of differences

plot(sigmaDIFF$sigmaT, main="Differences of sigma_t", 
     sub="sigmaLMER - sigmaSTAN",
     ylab="Difference of sigma_t Estimates")
abline(h=0, col="red")

plot(sigmaDIFF$sigmaP, main="Differences of sigma_p", 
     sub="sigmaLMER - sigmaSTAN",
     ylab="Difference of sigma_p Estimates")
abline(h=0, col="red")

plot(sigmaDIFF$sigmaE, main="Differences of sigma_e", 
     sub="sigmaLMER - sigmaSTAN",
     ylab="Difference of sigma_e Estimates")
abline(h=0, col="red")


# Plot of difference of var sum
plot(sigmaDIFF$var_sum, main="Var Sum of lmer - Var Sum of Stan", 
     ylab="Difference of Sum")
abline(h=0, col="red")

# Distribution of sigma
par(mfrow=c(1,2))
hist(sigmaLMER_200$sigmaT, main="Distribution of sigmaT (lmer)",
     xlab="sigmaT (lmer)", breaks=200,xlim=c(0,2.5))
hist(sigmaSTAN_200$sigmaT, main="Distribution of sigmaT (Stan)",
     xlab="sigmaT (Stan)", breaks=300,xlim=c(0,2.5))

par(mfrow=c(1,2))
hist(sigmaLMER_200$sigmaP, main="Distribution of sigmaP (lmer)",
     xlab="sigmaP (lmer)", breaks=200,xlim=c(0,2))
hist(sigmaSTAN_200$sigmaP, main="Distribution of sigmaP (Stan)",
     xlab="sigmaP (Stan)", breaks=300,xlim=c(0,2))


par(mfrow=c(1,2))
hist(sigmaLMER_200$sigmaE, main="Distribution of sigmaE (lmer)",
     xlab="sigmaE (lmer)", breaks=200,xlim=c(0,1.2))
hist(sigmaSTAN_200$sigmaE, main="Distribution of sigmaE (Stan)",
     xlab="sigmaE (Stan)", breaks=300,xlim=c(0,1.2))
