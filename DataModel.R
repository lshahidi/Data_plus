# this is the first code for modelling our data with lme4
# will incorporate RStan model as well


### INITIALIZE

# install.packages("coda")
# install.packages("gtools")

library(lme4)
library(arm)
library(rstan)
library(coda)
library(gtools)

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

Data1 <- site(1)
Data2 <- site(2)
Data3 <- site(3)
Data4 <- site(4)
Data5 <- site(5)
Data6 <- site(6)
Data7 <- site(7)
Data8 <- site(8)
Data9 <- site(9)


### FIT WITH LMER

fit1 <- lmer(X1 ~ tInd + (tInd|patient), Data1)
fit2 <- lmer(X2 ~ tInd + (tInd|patient), Data2)
fit3 <- lmer(X3 ~ tInd + (tInd|patient), Data3)
fit4 <- lmer(X4 ~ tInd + (tInd|patient), Data4)
fit5 <- lmer(X5 ~ tInd + (tInd|patient), Data5)
fit6 <- lmer(X6 ~ tInd + (tInd|patient), Data6)
fit7 <- lmer(X7 ~ tInd + (tInd|patient), Data7)
fit8 <- lmer(X8 ~ tInd + (tInd|patient), Data8)
fit9 <- lmer(X9 ~ tInd + (tInd|patient), Data9)

# fitGland <- lmer(X1 ~ tInd + (tInd|patient) + (0+gInd|patient:side), Data1)

patLabs <- c("C", "D","E","F","H","J","K","K*","M","O","P","S","T","U","W","X")

barplot(fixef(fit1), main="Fixed Effects", xlab="Effect", ylab="Intercept Estimate", names.arg=c("mu","betaT"), ylim=c(-2,2))
barplot(ranef(fit1)$patient$'(Intercept)', main="Patient Random Effect", xlab="Patient", ylab=expression(paste("Intercept Estimate (",b[i],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit1)$patient$tInd, main="Patient,Tumor Slope Random Effect", xlab="Patient", ylab=expression(paste("Slope Estimate (",b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit1)$patient$'(Intercept)'+ranef(fit1)$patient$tInd, main="Both Random Effects", xlab="Patient", ylab=expression(paste("Both Estimates (",b[i]+b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))

barplot(fixef(fit4), main="Fixed Effects", xlab="Effect", ylab="Intercept Estimate", names.arg=c("mu","betaT"), ylim=c(-2,2))
barplot(ranef(fit4)$patient$'(Intercept)', main="Patient Random Effect", xlab="Patient", ylab=expression(paste("Intercept Estimate (",b[i],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit4)$patient$tInd, main="Patient,Tumor Slope Random Effect", xlab="Patient", ylab=expression(paste("Slope Estimate (",b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit4)$patient$'(Intercept)'+ranef(fit1)$patient$tInd, main="Both Random Effects", xlab="Patient", ylab=expression(paste("Both Estimates (",b[i]+b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))


## extract variance at 1000 sites
nSites<-1000
variances <- data.frame(sigmaT=numeric(nSites),sigmaP=numeric(nSites),sigmaE=numeric(nSites))

ptm <- proc.time()
for (i in 1:nSites){
  print(paste("site ",i))
  
  DataI <- site(i)
  names(DataI)[1] <- "beta"
  
  fit <- lmer(beta ~ tInd + (tInd|patient), DataI)
  
  variances$sigmaT[i] <- as.data.frame(VarCorr(fit))$sdcor[1]
  variances$sigmaP[i] <- as.data.frame(VarCorr(fit))$sdcor[2]
  variances$sigmaE[i] <- as.data.frame(VarCorr(fit))$sdcor[4]
}
proc.time() - ptm

PTratio <- variances$sigmaP/variances$sigmaT
PTratio[PTratio>100] <- NA
hist(PTratio,100)

logPTratio <- log(variances$sigmaP/variances$sigmaT)
logPTratio[logPTratio>100] <- NA
hist(logPTratio,100)

library(reshape2)
library(ggplot2)
gg <- melt(variances)
ggplot(gg, aes(x=value, fill=variable)) +
  geom_histogram(binwidth=0.05)+
  facet_grid(variable~.)



# Extract standard deviation in lmer models

sd <- as.data.frame(VarCorr(fit1))[-3,5]
names(sd) <- c("sigma_p", "sigma_t", "sigma_e")
sd



# ALL SITES
# Test on first ten sites

Data <- list(NULL)
for (i in 1:10) {
  Data[[i]] <- site(i)
}

fit <- list(NULL)
for (i in 1:10) {
  fit[[i]] <- lmer(Data[[i]][,1] ~ tInd + (tInd|patient), Data[[i]])
}


temp <- as.data.frame(sapply(lapply(lapply(fit, VarCorr),as.data.frame), '[', i=5))[-3,]

sd_lmer <- t(temp)
colnames(sd_lmer) <- c("sigma_p", "sigma_t", "sigma_e")
rownames(sd_lmer) <- paste("Site",1:10)

sd_lmer



### FIT WITH STAN


# Function to build stan model for each site
stanfit <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])
  
  
  # Model 1: uniform priors
   stanFit <- stan(file="model.stan",data=stanDat)
  
  
  # Using Model 2: wide range uniform prior on sd and flat normal on fixed effects
  stanFit2 <- stan(file="model2.stan", data=stanDat)
  
  # Estimates for first three samples
   est <- head(as.matrix(stanFit,pars=c("mu","betaT","b_pat","bT_pat"))[,c(1:5,19:21)])

  est2 <- head(as.matrix(stanFit2,pars=c("mu","betaT","b_pat","bT_pat"))[,c(1:5,19:21)])

  # Variance
   var <- head(as.matrix(stanFit,pars=c("sigma_e","sigma_p","sigma_t")))
 
  var2 <- head(as.matrix(stanFit2,pars=c("sigma_e","sigma_p","sigma_t")))
 
  return(list(stanFit=stanFit, stanFit2=stanFit2, est=est, est2=est2, var=var, var2=var2))
}

stan1 <- stanfit(Data1)
summary(stan1$stanFit2)

stan2 <- stanfit(Data2)
summary(stan2$stanFit2)


stan3 <- stanfit(Data3)
summary(stan3$stanFit2)

stan4 <- stanfit(Data4)
summary(stan4$stanFit2)

stan5 <- stanfit(Data5)
summary(stan5$stanFit2)

stan6 <- stanfit(Data6)
summary(stan6$stanFit2)

stan7 <- stanfit(Data7)
summary(stan7$stanFit2)

stan8 <- stanfit(Data8)
summary(stan8$stanFit2)


plot(stan1$stanFit2, pars=c("sigma_t","sigma_p","sigma_e"))
plot(stan2$stanFit2, pars=c("sigma_t","sigma_p","sigma_e"))
plot(stan3$stanFit2, pars=c("sigma_t","sigma_p","sigma_e"))

# site 4: sigma_t < sigma_p
plot(stan4$stanFit2, pars=c("sigma_t","sigma_p","sigma_e"))

plot(stan5$stanFit2, pars=c("sigma_t","sigma_p","sigma_e"))
plot(stan6$stanFit2, pars=c("sigma_t","sigma_p","sigma_e"))
plot(stan7$stanFit2, pars=c("sigma_t","sigma_p","sigma_e"))
plot(stan8$stanFit2, pars=c("sigma_t","sigma_p","sigma_e"))




# Extract standard deviation from stan models
summary(stan1$stanFit2)$summary[35:37,1]



# ALL SITES
# Test on first ten sites

stan <- lapply(Data, stanfit)

for (i in 1:10) {
  rbind(sd_stan,summary(stan[i]$stanFit2)$summary[35:37,1])
}




# Diagnostic Graphs
mcmcCoda <- mcmc.list(lapply( 1:ncol(stanFit) ,
                              function(x) { mcmc(as.array(stanFit)[,x,])} ))

plot(stanFit, pars=c("sigma_t","sigma_p","sigma_e"), main="Variances (Model 1)")
plot(stanFit2, pars=c("sigma_t","sigma_p","sigma_e"), main="Variances (Model 2)")



