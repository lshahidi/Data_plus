############################ Stan Model Comparison ##############################

library(rstan)
library(gtools)
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


load("betaT.RData")
load("betaT2.RData")
load("mu.RData")
load("mu2.RData")
load("sigma_t.RData")
load("sigma_t2.RData")
load("sigma_p.RData")
load("sigma_p1.RData")
load("sigma_pt2.RData")
load("sigma_e.RData")
load("sigma_e2.RData")


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

# We have full results from model 3
# To solve all the warnings, tried model_rep, model_rep2, model_try, model_try2
# Only model_rep return no divergence warnings


# Using Model 2: wide range uniform prior on sd and flat normal on fixed effects

stanfit2 <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])
  
  
  
  stanFit2 <- stan(file="model2.stan", data=stanDat, control = list(adapt_delta = 0.999))
  
  
  return(stanFit2=stanFit2)
}


# Using Model 3: add intra-tumoral variance

stanfit3 <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])
  
  
  
  stanFit3 <- stan(file="model3.stan", data=stanDat, control = list(adapt_delta = 0.999))
  
  return(stanFit3=stanFit3)
}


# Using Model 4: Simple Model - only use tumor samples

stanfit4 <- function (dataset) {
  
  data <- dataset[dataset[,3]=="T",]
  
  stanDat <- list(pID = as.integer(factor(data$patient)),
                  N = nrow(data),
                  P = nlevels(data$patient),
                  y = data[,1])
  
  
  stanFit4 <- stan(file="model4.stan", data=stanDat, control = list(adapt_delta = 0.999))
  
  return(stanFit4=stanFit4)
}


# Using Reparameterized Model (Fully reparameterized)

stanfit_rep <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])
  
  
  # Using reparameterized model
  
  stanFit_rep <- stan(file="model_rep.stan", data=stanDat, control=list(adapt_delta=0.999))
  
  return(stanFit=stanFit_rep)
}


# Using reparameterized model (reparameterization only on b_pat, c_patT, d_T)

stanfit_rep2 <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])
  
  
  
  stanFit_rep2 <- stan(file="model_rep2.stan", data=stanDat, control=list(adapt_delta=0.999))
  
  return(stanFit=stanFit_rep2)
}


# Try: gamma(2,2) prior on sigma, partially reparameterized

stanfit_try <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])
  
  
  
  stanFit_try <- stan(file="model_try.stan", data=stanDat, control=list(adapt_delta=0.999))
  
  return(stanFit=stanFit_try)
}


# Try2: gamma(2,2) prior on sigma, fully reparameterized

stanfit_try2 <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])
  
  
  
  stanFit_try2 <- stan(file="model_try2.stan", data=stanDat, control=list(adapt_delta=0.999))
  
  return(stanFit=stanFit_try2)
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
Data10 <- site(10)




stan1 <- stanfit3(Data1)
stan1_rep <- stanfit_rep(Data1)
stan1_rep2 <- stanfit_rep2(Data1)
stan1_try <- stanfit_try(Data1)


stan2 <- stanfit3(Data2)
stan2_rep <- stanfit_rep(Data2)
stan2_rep2 <- stanfit_rep2(Data2)
stan2_try <- stanfit_try(Data2)



stan3 <- stanfit3(Data3)
stan3_rep <- stanfit_rep(Data3)
stan3_rep2 <- stanfit_rep2(Data3)
stan3_try <- stanfit_try(Data3)


stan4 <- stanfit3(Data4)
stan4_rep <- stanfit_rep(Data4)
stan4_rep2 <- stanfit_rep2(Data4)
stan4_try <- stanfit_try(Data4)




stan5 <- stanfit3(Data5)
stan5_rep <- stanfit_rep(Data5)
stan5_rep2 <- stanfit_rep2(Data5)
stan5_try <- stanfit_try(Data5)



stan6 <- stanfit3(Data6)
stan6_rep <- stanfit_rep(Data6)
stan6_rep2 <- stanfit_rep2(Data6)
stan6_try <- stanfit_try(Data6)


stan7 <- stanfit3(Data7)
stan7_rep <- stanfit_rep(Data7)
stan7_rep2 <- stanfit_rep2(Data7)
stan7_try <- stanfit_try(Data7)



Data39490 <- site(39490)
stan39490_rep <- stanfit_rep(Data39490)


test <- Data2

temp <- data.frame(X2=c(-1.44,2,1.24,-0.994,-0.83,0.66,1.24,-0.994,-0.8,1.02,-0.98,0.66),
                   patient=c("A","A","A","B","B","B","R","R","R","Y","Y","Y"),
                   tissue=c("T","T","N","T","T","N","T","T","N","T","T","N"),
                   side=c("A","B","N","A","B","N","A","B","N","A","B","N"),
                   tInd=c(1,1,0,1,1,0,1,1,0,1,1,0))
test <- rbind(test,temp)

fit <- stanfit3(test)

stan1 <- stanfit3(Data1)
stan2 <- stanfit3(Data2)
stan3 <- stanfit4(Data3)
stan4 <- stanfit4(Data4)

stan1 <- stanfit2(Data1)
stan1t <- stanfit3(Data1)

stan2 <- stanfit2(Data2)
stan2t <- stanfit3(Data2)

stan3 <- stanfit2(Data3)
stan3t <- stanfit3(Data3)

stan4 <- stanfit2(Data4)
stan4t <- stanfit3(Data4)

plot(stan1, pars=c("sigma_t","sigma_p","sigma_e"))
plot(stan1t, pars=c("sigma_t","sigma_p","sigma_e","sigma_pt"))

plot(stan2, pars=c("sigma_t","sigma_p","sigma_e"))
plot(stan2t, pars=c("sigma_t","sigma_p","sigma_e","sigma_pt"))

plot(stan3, pars=c("sigma_t","sigma_p","sigma_e"))
plot(stan3t, pars=c("sigma_t","sigma_p","sigma_e","sigma_pt"))

plot(stan4, pars=c("sigma_t","sigma_p","sigma_e"))
plot(stan4t, pars=c("sigma_t","sigma_p","sigma_e","sigma_pt"))



# Explore posterior distributions

posterior1 <- as.array(stan1)

mcmc_areas(
  posterior1, 
  pars = c("sigma_p", "sigma_pt","sigma_t","sigma_e"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, 
  point_est = "median"
)

posterior1t <- as.array(stan1t)

mcmc_areas(
  posterior1t, 
  pars = c("sigma_e", "sigma_p", "sigma_t","sigma_pt"),
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

posterior2t <- as.array(stan2t)

mcmc_areas(
  posterior2t, 
  pars = c("sigma_e", "sigma_p", "sigma_t","sigma_pt"),
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

posterior3t <- as.array(stan3t)

mcmc_areas(
  posterior3t, 
  pars = c("sigma_e", "sigma_p", "sigma_t","sigma_pt"),
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

posterior4t <- as.array(stan4t)

mcmc_areas(
  posterior4t, 
  pars = c("sigma_e", "sigma_p", "sigma_t","sigma_pt"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, 
  point_est = "mean"
)




# Posterior ratio and probability

prob <- rep(NA,dim(FullAnnotation)[1])

for (i in dim(FullAnnotation)[1]) {
  
  data <- site(i)
  stan <- stanfit3(data)

  posterior <- as.matrix(stan,pars=c("sigma_p","sigma_t"))

  PTratio <- log(sigma[,1]/sigma[,2])

  prob[i] <- sum(PTratio1 > 0)/length(PTratio)

}






# Extract estimates from stan models on random 200 sites

set.seed(123)
indice <- sample(1:dim(FullAnnotation)[1],200,replace=FALSE)
indice <- sort(indice)
nsites <- 200


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


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

count <- 1

for (i in indice) {
  data <- site(i)
  stan <- stanfit2(data)
  posterior <- as.array(stan)
  
  betaT[count,1] <- summary(stan)$summary[33,1]
  betaT[count,2] <- getmode(posterior[1:1000,1:4,33])
  betaT[count,3:7] <- summary(stan)$summary[33,4:8]
  
  mu[count,1] <- summary(stan)$summary[34,1]
  mu[count,2] <- getmode(posterior[1:1000,1:4,34])
  mu[count,3:7] <- summary(stan)$summary[34,4:8]
  
  sigma_e[count,1] <- summary(stan)$summary[35,1]
  sigma_e[count,2] <- getmode(posterior[1:1000,1:4,35])
  sigma_e[count,3:7] <- summary(stan)$summary[35,4:8]
  
  
  sigma_p[count,1] <- summary(stan)$summary[36,1]
  sigma_p[count,2] <- getmode(posterior1[1:1000,1:4,36])
  sigma_p[count,3:7] <- summary(stan)$summary[36,4:8]
  
  
  sigma_t[count,1] <- summary(stan)$summary[37,1]
  sigma_t[count,2] <- getmode(posterior1[1:1000,1:4,37])
  sigma_t[count,3:7] <- summary(stan)$summary[37,4:8]
  
  count <- count + 1
}




betaT2 <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                    p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                    p97.5=numeric(nsites))

mu2 <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                 p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                 p97.5=numeric(nsites))

sigma_e2 <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                      p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                      p97.5=numeric(nsites))


sigma_p2 <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                      p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                      p97.5=numeric(nsites))


sigma_t2 <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                      p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                      p97.5=numeric(nsites))


sigma_pt2 <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                       p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                       p97.5=numeric(nsites))

count <- 1

for (i in indice) {
  data <- site(i)
  stan <- stanfit3(data)
  posterior <- as.array(stan)
  
  betaT2[count,1] <- summary(stan)$summary[71,1]
  betaT2[count,2] <- getmode(posterior[1:1000,1:4,71])
  betaT2[count,3:7] <- summary(stan)$summary[71,4:8]
  
  mu2[count,1] <- summary(stan)$summary[72,1]
  mu2[count,2] <- getmode(posterior[1:1000,1:4,72])
  mu2[count,3:7] <- summary(stan)$summary[72,4:8]
  
  sigma_e2[count,1] <- summary(stan)$summary[73,1]
  sigma_e2[count,2] <- getmode(posterior[1:1000,1:4,73])
  sigma_e2[count,3:7] <- summary(stan)$summary[73,4:8]
  
  
  sigma_p2[count,1] <- summary(stan)$summary[74,1]
  sigma_p2[count,2] <- getmode(posterior[1:1000,1:4,74])
  sigma_p2[count,3:7] <- summary(stan)$summary[74,4:8]
  
  sigma_pt2[count,1] <- summary(stan)$summary[75,1]
  sigma_pt2[count,2] <- getmode(posterior[1:1000,1:4,75])
  sigma_pt2[count,3:7] <- summary(stan)$summary[75,4:8]
  
  sigma_t2[count,1] <- summary(stan)$summary[76,1]
  sigma_t2[count,2] <- getmode(posterior[1:1000,1:4,76])
  sigma_t2[count,3:7] <- summary(stan)$summary[76,4:8]
  
  count <- count + 1
}



save(betaT, file="betaT.RData")
save(betaT2,file="betaT2.RData")
save(mu, file="mu.RData")
save(mu2,file="mu2.RData")
save(sigma_t, file="sigma_t.RData")
save(sigma_t2, file="sigma_t2.RData")
save(sigma_p, file="sigma_p.RData")
save(sigma_p2, file="sigma_p2.RData")
save(sigma_e, file="sigma_e.RData")
save(sigma_e2, file="sigma_e2.RData")
save(sigma_pt2, file="sigma_pt2.RData")



hist(mu$mean, main="Distribution of mu",
     breaks = 100,xlim=c(-5,5),col=grey.colors(1,alpha=0.1))
hist(mu2$mean,breaks=100,col="black",add=TRUE)
legend("topright",c("New","Old"),fill=c(grey.colors(1,alpha=0.1),"black"))


hist(betaT$mean, main="Distribution of betaT",
     breaks = 100,xlim=c(-3,3),col=grey.colors(1,alpha=0.1))
hist(betaT2$mean,breaks=100,col="black",add=TRUE)
legend("topright",c("New","Old"),fill=c(grey.colors(1,alpha=0.1),"black"))


hist(sigma_p$mean, main="Distribution of sigmaP",
     breaks = 100,xlim=c(0,1.5),col=grey.colors(1,alpha=0.1))
hist(sigma_p2$mean,breaks=100,col="black",add=TRUE)
legend("topright",c("New","Old"),fill=c(grey.colors(1,alpha=0.1),"black"))


hist(sigma_t$mean, main="Distribution of sigmaT (at patient level)",
     breaks = 100,xlim=c(0,2.5),col=grey.colors(1,alpha=0.1))
hist(sigma_pt2$mean,breaks=100,col="black",add=TRUE)
legend("topright",c("New","Old"),fill=c(grey.colors(1,alpha=0.1),"black"))


hist(sigma_t2$mean, main="Distribution of Intratumoral SD",
     breaks = 100,xlim=c(0,0.2),col="black")


hist(sigma_e$mean, main="Distribution of sigmaE",
     breaks = 100,xlim=c(0,1.5),col=grey.colors(1,alpha=0.1))
hist(sigma_e2$mean,breaks=100,col="black",add=TRUE)
legend("topright",c("New","Old"),fill=c(grey.colors(1,alpha=0.1),"black"))





####### Compare Complex and Simple Stan Models Using 5K Sites #########

set.seed(555)
indice <- sample(1:dim(FullAnnotation)[1],5000,replace=FALSE)
indice <- sort(indice)
nsites <- 5000


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


betaT_C <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                     p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                     p97.5=numeric(nsites))

mu_C <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                  p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                  p97.5=numeric(nsites))

sigmaE_C <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                       p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                       p97.5=numeric(nsites))


sigmaP_C <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                       p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                       p97.5=numeric(nsites))


sigmaT_C <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                       p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                       p97.5=numeric(nsites))


sigmaPT_C <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                        p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                        p97.5=numeric(nsites))

count <- 1

for (i in indice) {
  data <- site(i)
  stan <- stanfit3(data)
  posterior <- as.array(stan)
  
  betaT_C[count,1] <- summary(stan)$summary[71,1]
  betaT_C[count,2] <- getmode(posterior[1:1000,1:4,71])
  betaT_C[count,3:7] <- summary(stan)$summary[71,4:8]
  
  mu_C[count,1] <- summary(stan)$summary[72,1]
  mu_C[count,2] <- getmode(posterior[1:1000,1:4,72])
  mu_C[count,3:7] <- summary(stan)$summary[72,4:8]
  
  sigmaE_C[count,1] <- summary(stan)$summary[73,1]
  sigmaE_C[count,2] <- getmode(posterior[1:1000,1:4,73])
  sigmaE_C[count,3:7] <- summary(stan)$summary[73,4:8]
  
  
  sigmaP_C[count,1] <- summary(stan)$summary[74,1]
  sigmaP_C[count,2] <- getmode(posterior[1:1000,1:4,74])
  sigmaP_C[count,3:7] <- summary(stan)$summary[74,4:8]
  
  sigmaPT_C[count,1] <- summary(stan)$summary[75,1]
  sigmaPT_C[count,2] <- getmode(posterior[1:1000,1:4,75])
  sigmaPT_C[count,3:7] <- summary(stan)$summary[75,4:8]
  
  sigmaT_C[count,1] <- summary(stan)$summary[76,1]
  sigmaT_C[count,2] <- getmode(posterior[1:1000,1:4,76])
  sigmaT_C[count,3:7] <- summary(stan)$summary[76,4:8]
  
  count <- count + 1
}



# Partial Results
load("mu_C.RData")
load("mu_S.RData")
load("betaT_C.RData")
load("sigmaE_C.RData")
load("sigmaE_S.RData")
load("sigmaP_C.RData")
load("sigmaP_S.RData")
load("sigmaT_C.RData")
load("sigmaT_S.RData")
load("sigmaPT_C.RData")

hist(mu_C[1:621,1],main="mu",xlab="mu",
     breaks=300,col=grey.colors(1,alpha=0.1))
hist(mu_S[1:621,1],breaks=300,col="coral",add=TRUE)
legend("topright",c("Complex","Simple"),fill=c(grey.colors(1,alpha=0.1),"coral"))

hist(betaT_C[1:621,1],main="betaT (Only for Complex Model)",xlab="betaT",
     breaks=300,col=grey.colors(1,alpha=0.1))

hist(sigmaP_C[1:621,1],main="sigmaP",xlab="sigmaP",
     breaks=300,col=grey.colors(1,alpha=0.1))
hist(sigmaP_S[1:621,1],breaks=300,col="coral",add=TRUE)
legend("topright",c("Complex","Simple"),fill=c(grey.colors(1,alpha=0.1),"coral"))

hist(sigmaPT_C[1:621,1],main="sigmaPT (Only for Complex Model)",xlab="sigmaPT",
     breaks=300,col=grey.colors(1,alpha=0.1))

hist(sigmaT_C[1:621,1],main="sigmaT",xlab="sigmaT",
     breaks=300,col=grey.colors(1,alpha=0.1))
hist(sigmaT_S[1:621,1],breaks=300,col="coral",add=TRUE)
legend("topright",c("Complex","Simple"),fill=c(grey.colors(1,alpha=0.1),"coral"))


hist(sigmaE_C[1:621,1],main="sigmaE",xlab="sigmaE",
     breaks=300,col=grey.colors(1,alpha=0.1))
hist(sigmaE_S[1:621,1],breaks=300,col=adjustcolor("coral",0.5),
     add=TRUE)
legend("topright",c("Complex","Simple"),fill=c(grey.colors(1,alpha=0.1),adjustcolor("coral",0.5)))



hist((log((sigmaP_C$mean/sigmaT_C$mean)))[1:621],main="logPTratio",
     xlab="log sigmaP/sigmaT", xlim=c(-2.5,3.5), ylim=c(0,20),
     breaks=300, col=grey.colors(1,alpha=0.1))
hist((log((sigmaP_S$mean/sigmaT_S$mean)))[1:621],breaks=300,
     col=adjustcolor("black",0.5),add=TRUE)
legend("topright",c("Complex","Simple"),fill=c(grey.colors(1,alpha=0.1),"black"))



######### Compare original model with reparameterized model (random 200 sites) ############

set.seed(123)
indice <- sample(1:dim(FullAnnotation)[1],200,replace=FALSE)
indice <- sort(indice)
nsites <- 200

# save(indice,file="indice.RData")

load("indice.RData")



est <- data.frame(betaT=numeric(nsites),mu=numeric(nsites),sigma_e=numeric(nsites),
                    sigma_p=numeric(nsites),sigma_pt=numeric(nsites),sigma_t=numeric(nsites))


count <- 1

for (i in indice) {
  data <- site(i)
  stan <- stanfit3(data)
  
  est[count,] <- summary(stan)$summary[71:76,1]

  count <- count + 1
}



est_rep3 <- data.frame(betaT=numeric(nsites),mu=numeric(nsites),sigma_e=numeric(nsites),
                  sigma_p=numeric(nsites),sigma_pt=numeric(nsites),sigma_t=numeric(nsites))

rownames(est_rep3) <- paste("site",indice)

count <- 1

for (i in indice) {
  data <- site(i)
  stan <- stanfit_rep(data)
  
  est_rep3[count,] <- summary(stan)$summary[147:152,1]
  
  count <- count + 1
}

load("rep.RData")
load("norep.RData")

DIFF <- est - est_rep

plot(DIFF$mu, main="Difference of mu", 
     ylab="Original Model - Reparameterized Model",ylim=c(-15,15))
abline(h=0,col="red")
text(175,10,labels=paste("Mean = ",round(mean(DIFF$mu),4)))

plot(DIFF$betaT, main="Difference of betaT", 
     ylab="Original Model - Reparameterized Model", ylim=c(-15,15))
abline(h=0,col="red")
text(175,10,labels=paste("Mean = ",round(mean(DIFF$betaT),4)))

plot(DIFF$sigma_p, main="Difference of sigmaP", 
     ylab="Original Model - Reparameterized Model")
abline(h=0,col="red")
text(175,-15,labels=paste("Mean = ",round(mean(DIFF$sigma_p),4)))


plot(DIFF$sigma_pt, main="Difference of sigmaPT", 
     ylab="Original Model - Reparameterized Model")
abline(h=0,col="red")
text(175,-15,labels=paste("Mean = ",round(mean(DIFF$sigma_pt),4)))


plot(DIFF$sigma_t, main="Difference of sigmaT", 
     ylab="Original Model - Reparameterized Model")
abline(h=0,col="red")
text(175,-0.5,labels=paste("Mean = ",round(mean(DIFF$sigma_t),4)))

plot(DIFF$sigma_e, main="Difference of sigmaE", 
     ylab="Original Model - Reparameterized Model", ylim=c(-1.2,1))
abline(h=0,col="red")
text(175,0.5,labels=paste("Mean = ",round(mean(DIFF$sigma_e),4)))



#### Check the Variablity of Estimates from Reparameterized Model 

load("rep.RData")
load("rep2.RData")
load("rep3.RData")


DIFF12 <- est_rep - est_rep2
DIFF13 <- est_rep - est_rep3
DIFF23 <- est_rep2 - est_rep3

var12 <- DIFF12/est_rep2
var13 <- DIFF13/est_rep3
var23 <- DIFF23/est_rep3

plot(var12$mu,main="Variablity of mu",ylab="Change of mu")
text(25,-200,labels=paste("Mean = ",round(mean(var12$mu),2)))
plot(var12$mu,main="Variablity of mu",ylab="Change of mu",ylim=c(-10,15))
abline(h=c(-0.2,0.2),col="red")
text(25,13,labels=paste(">20% = ", sum(abs(var12$mu)>0.2)))


plot(var12$betaT,main="Variablity of betaT",ylab="Change of betaT")
text(25,1000,labels=paste("Mean = ",round(mean(var12$betaT),2)))
plot(var12$mu,main="Variablity of mu",ylab="Change of mu",ylim=c(-10,15))
abline(h=c(-0.2,0.2),col="red")
text(25,13,labels=paste(">20% = ", sum(abs(var12$betaT)>0.2)))


plot(var12$sigma_p,main="Variablity of sigmaP",ylab="Change of sigmaP")
text(100,35,labels=paste("Mean = ",round(mean(var12$sigma_p),2)))
abline(h=c(-0.2,0.2),col="red")
text(100,30,labels=paste(">20% = ", sum(abs(var12$sigma_p)>0.2)))


plot(var12$sigma_pt,main="Variablity of sigmaPT",ylab="Change of sigmaPT")
text(100,45,labels=paste("Mean = ",round(mean(var12$sigma_pt),2)))
abline(h=c(-0.2,0.2),col="red")
text(100,40,labels=paste(">20% = ", sum(abs(var12$sigma_pt)>0.2)))


plot(var12$sigma_t,main="Variablity of sigmaT",ylab="Change of sigmaT")
text(100,2.5,labels=paste("Mean = ",round(mean(var12$sigma_t),2)))
abline(h=c(-0.2,0.2),col="red")
text(100,2,labels=paste(">20% = ", sum(abs(var12$sigma_t)>0.2)))


plot(var12$sigma_e,main="Variablity of sigmaE",ylab="Change of sigmaE")
text(25,0.6,labels=paste("Mean = ",round(mean(var12$sigma_e),2)))
abline(h=c(-0.2,0.2),col="red")
text(25,0.5,labels=paste(">20% = ", sum(abs(var12$sigma_e)>0.2)))


var12$varno <- apply(abs(var12)>0.2,1,sum)
var12$flag <- ifelse(var12$varno>0,1,0)
var23$varno <- apply(abs(var23)>0.2,1,sum)
var23$flag <- ifelse(var23$varno>0,1,0)
var13$varno <- apply(abs(var13)>0.2,1,sum)
var13$flag <- ifelse(var13$varno>0,1,0)

flag <- rep(NA,200)
for (i in 1:200) {
  if (var12$flag[i] == 1 | var23$flag[i] == 1 | var13$flag[i] == 1) {
    flag[i] = 1
  }
}




#### Check the Variablity of Estimates from Non-reparameterized Model 

load("norep.RData")
load("norep2.RData")
load("norep3.RData")


DIFF12 <- est - est2
DIFF13 <- est - est3
DIFF23 <- est2 - est3

var12 <- DIFF12/est2
var13 <- DIFF13/est3
var23 <- DIFF23/est3


plot(var12$mu,main="Variablity of mu",ylab="Change of mu")
abline(h=c(-0.2,0.2),col="red")
text(175,0.7,labels=paste(">20% = ", sum(abs(var12$mu)>0.2)))
abline(h=c(-0.05,0.05),col="blue")
text(175,0.6,labels=paste(">5% = ", sum(abs(var12$mu)>0.05)))
legend(160,0.5,lty=c(1,1),col=c("red","blue"),legend=c("20%","5"),cex=0.8)


plot(var12$betaT,main="Variablity of betaT",ylab="Change of betaT")
plot(var12$betaT,main="Variablity of betaT",ylab="Change of betaT",ylim=c(-1.5,1))
abline(h=c(-0.2,0.2),col="red")
text(175,-0.6,labels=paste(">20% = ", sum(abs(var12$betaT)>0.2)))
abline(h=c(-0.05,0.05),col="blue")
text(175,-0.8,labels=paste(">5% = ", sum(abs(var12$betaT)>0.05)))
legend(160,-1,lty=c(1,1),col=c("red","blue"),legend=c("20%","5"),cex=0.7)


plot(var12$sigma_p,main="Variablity of sigmaP",ylab="Change of sigmaP")
abline(h=c(-0.2,0.2),col="red")
text(50,0.6,labels=paste(">20% = ", sum(abs(var12$sigma_p)>0.2)))
abline(h=c(-0.05,0.05),col="blue")
text(50,0.45,labels=paste(">5% = ", sum(abs(var12$sigma_p)>0.05)))
legend(180,0.5,lty=c(1,1),col=c("red","blue"),legend=c("20%","5"),cex=0.7)



plot(var12$sigma_pt,main="Variablity of sigmaPT",ylab="Change of sigmaPT",ylim=c(-0.3,0.3))
abline(h=c(-0.2,0.2),col="red")
text(25,0.28,labels=paste(">20% = ", sum(abs(var12$sigma_pt)>0.2)))
abline(h=c(-0.05,0.05),col="blue")
text(60,0.28,labels=paste(">5% = ", sum(abs(var12$sigma_pt)>0.05)))
legend(180,0.3,lty=c(1,1),col=c("red","blue"),legend=c("20%","5"),cex=0.7)



plot(var12$sigma_t,main="Variablity of sigmaT",ylab="Change of sigmaT")
abline(h=c(-0.2,0.2),col="red")
text(175,0.9,labels=paste(">20% = ", sum(abs(var12$sigma_t)>0.2)))
abline(h=c(-0.05,0.05),col="blue")
text(175,0.75,labels=paste(">5% = ", sum(abs(var12$sigma_t)>0.05)))
legend(170,0.65,lty=c(1,1),col=c("red","blue"),legend=c("20%","5"),cex=0.7)





plot(var12$sigma_e,main="Variablity of sigmaE",ylab="Change of sigmaE")
abline(h=c(-0.2,0.2),col="red")
text(25,0.65,labels=paste(">20% = ", sum(abs(var12$sigma_e)>0.2)))
abline(h=c(-0.05,0.05),col="blue")
text(25,0.5,labels=paste(">5% = ", sum(abs(var12$sigma_e)>0.05)))
legend(180,0.65,lty=c(1,1),col=c("red","blue"),legend=c("20%","5"),cex=0.7)




var12$varno <- apply(abs(var12)>0.2,1,sum)
var12$flag <- ifelse(var12$varno>0,1,0)
var23$varno <- apply(abs(var23)>0.2,1,sum)
var23$flag <- ifelse(var23$varno>0,1,0)
var13$varno <- apply(abs(var13)>0.2,1,sum)
var13$flag <- ifelse(var13$varno>0,1,0)

flag <- rep(NA,200)
for (i in 1:200) {
  if (var12$flag[i] == 1 | var23$flag[i] == 1 | var13$flag[i] == 1) {
    flag[i] = 1
  } else {
    flag[i] = 0
  }
}