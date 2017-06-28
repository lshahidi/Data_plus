# script for extracting variance from Stan model 
# part 1: on sites 1-425K


library(rstan)
library(gtools)

# load data
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


# Function to build stan model for each site
stanfit <- function (dataset) {
  
  stanDat <- list(pID = as.integer(factor(dataset$patient)),
                  tInd = dataset$tInd,
                  N = nrow(dataset),
                  P = nlevels(dataset$patient),
                  y = dataset[,1])

  # Using Model 2: wide range uniform prior on sd and flat normal on fixed effects
  stanFit2 <- stan(file="model2.stan", data=stanDat, control = list(adapt_delta = 0.999))

  return(stanFit2=stanFit2)
}

# Function to get mode

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


nsites <- dim(FullAnnotation)[1]/2

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


for (i in 1:nsites) {
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




# save data
save(betaT,file="betaTSTAN1.RData")
save(mu,file="muSTAN1.RData")
save(sigma_e,file="sigmaESTAN1.RData")
save(sigma_p,file="sigmaPSTAN1.RData")
save(sigma_t,file="sigmaTSTAN1.RData")

