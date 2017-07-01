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
stanfit2 <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])
  
  
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
  pars = c("sigma_e", "sigma_p", "sigma_t"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, 
  point_est = "mean"
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

