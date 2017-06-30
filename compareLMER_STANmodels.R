######################## Stan and lmer Comparison ################################



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




# Randomly choose 200 sites

set.seed(123)
indice <- sample(1:dim(FullAnnotation)[1],200,replace=FALSE)
indice <- sort(indice)


################### Old model (without intratumor variance) #####################

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



################### New model (with intratumor variance) #####################

# lmer
nSites<-200
sigmaLMER_new <- data.frame(sigmaT=numeric(nSites),sigmaTP=numeric(nSites),sigmaP=numeric(nSites),sigmaE=numeric(nSites),betaT=numeric(nSites),mu=numeric(nSites))

for (i in indice){
  print(paste("site ",i))
  
  DataI <- site(i)
  names(DataI)[1] <- "beta"
  
  fit <- lmer(beta ~ tInd + (1|patient) + (tInd|patient), DataI)
  
  sigmaLMER$sigmaT[i] <- as.data.frame(VarCorr(fit))$sdcor[1]
  sigmaLMER$sigmaP[i] <- as.data.frame(VarCorr(fit))$sdcor[2]
  sigmaLMER$sigmaE[i] <- as.data.frame(VarCorr(fit))$sdcor[4]
  sigmaLMER$betaT[i] <- fixef(fit)[2]
  sigmaLMER$mu[i] <- fixef(fit)[1]
}


