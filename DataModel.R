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

fitGland <- lmer(X1 ~ tInd + (tInd|patient) + (0+gInd|patient:side), Data1)

patLabs <- c("C", "D","E","F","H","J","K","K*","M","O","P","S","T","U","W","X")

barplot(fixef(fit1), main="Fixed Effects", xlab="Effect", ylab="Intercept Estimate", names.arg=c("mu","betaT"), ylim=c(-2,2))
barplot(ranef(fit1)$patient$'(Intercept)', main="Patient Random Effect", xlab="Patient", ylab=expression(paste("Intercept Estimate (",b[i],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit1)$patient$tInd, main="Patient,Tumor Slope Random Effect", xlab="Patient", ylab=expression(paste("Slope Estimate (",b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit1)$patient$'(Intercept)'+ranef(fit1)$patient$tInd, main="Both Random Effects", xlab="Patient", ylab=expression(paste("Both Estimates (",b[i]+b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))

barplot(fixef(fit2), main="Fixed Effects", xlab="Effect", ylab="Intercept Estimate", names.arg=c("mu","betaT"), ylim=c(-2,2))
barplot(fixef(fit22), main="Fixed Effects", xlab="Effect", ylab="Intercept Estimate", names.arg=c("mu"), ylim=c(-2,2))
barplot(ranef(fit2)$patient$'(Intercept)', main="Patient Random Effect", xlab="Patient", ylab=expression(paste("Intercept Estimate (",b[i],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit22)$patient$'(Intercept)', main="Patient Random Effect", xlab="Patient", ylab=expression(paste("Intercept Estimate (",b[i],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit2)$patient$tInd, main="Patient,Tumor Slope Random Effect", xlab="Patient", ylab=expression(paste("Slope Estimate (",b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit22)$patient$tInd, main="Patient,Tumor Slope Random Effect", xlab="Patient", ylab=expression(paste("Slope Estimate (",b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit2)$patient$'(Intercept)'+ranef(fit1)$patient$tInd, main="Both Random Effects", xlab="Patient", ylab=expression(paste("Both Estimates (",b[i]+b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit22)$patient$'(Intercept)'+ranef(fit1)$patient$tInd, main="Both Random Effects", xlab="Patient", ylab=expression(paste("Both Estimates (",b[i]+b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))


barplot(fixef(fit4), main="Fixed Effects", xlab="Effect", ylab="Intercept Estimate", names.arg=c("mu","betaT"), ylim=c(-2,2))
barplot(ranef(fit4)$patient$'(Intercept)', main="Patient Random Effect", xlab="Patient", ylab=expression(paste("Intercept Estimate (",b[i],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit4)$patient$tInd, main="Patient,Tumor Slope Random Effect", xlab="Patient", ylab=expression(paste("Slope Estimate (",b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))
barplot(ranef(fit4)$patient$'(Intercept)'+ranef(fit1)$patient$tInd, main="Both Random Effects", xlab="Patient", ylab=expression(paste("Both Estimates (",b[i]+b[iT],")")), names.arg = patLabs, cex.names = 0.7, ylim=c(-2,2))



### FIT WITH STAN


# Logit transform on beta values

#Data1$y <- log(Data1$X1/(1-Data1$X1))
Data1$y <- Data1$X1

# Stan data

stanDat <- list(pID = as.integer(factor(Data1$patient)),
                tInd = Data1$tInd,
                N = nrow(Data1),
                P = nlevels(Data1$patient),
                y = Data1$y)


stanFit <- stan(file="model.stan",data=stanDat)


# Using Model 2: wide range uniform prior on sd and flat normal on fixed effects
stanFit2 <- stan(file="model2.stan", data=stanDat)


print(stanFit, probs = c(0.025, 0.5, 0.975))

summary(stanFit)
summary(stanFit2)

est <- as.matrix(stanFit,pars=c("mu","betaT","b_pat","bT_pat"))
# Estimates for first three samples
head(est[,c(1:5,19:21)])


est2 <- as.matrix(stanFit2,pars=c("mu","betaT","b_pat","bT_pat"))
# Estimates for first three samples
head(est2[,c(1:5,19:21)])


# Variance
var <- as.matrix(stanFit,pars=c("sigma_e","sigma_p","sigma_t"))
head(var)

var2 <- as.matrix(stanFit2,pars=c("sigma_e","sigma_p","sigma_t"))
head(var2)




# Diagnostic Graphs
mcmcCoda <- mcmc.list(lapply( 1:ncol(stanFit) ,
                              function(x) { mcmc(as.array(stanFit)[,x,])} ))

par(mfrow=c(1,2))
plot(stanFit, pars=c("sigma_t","sigma_p","sigma_e"), main="Variances (Model 1)")
plot(stanFit2, pars=c("sigma_t","sigma_p","sigma_e"), main="Variances (Model 2)")










