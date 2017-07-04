# script for extracting variance from Stan model on random 5K sites
# SIMPLE MODEL
# part 6


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

# Function to build stan model for each site

stanfit4 <- function (dataset) {
  
  data <- dataset[dataset[,3]=="T",]
  
  stanDat <- list(pID = as.integer(factor(data$patient)),
                  N = nrow(data),
                  P = nlevels(data$patient),
                  y = data[,1])
  
  
  stanFit4 <- stan(file="model4.stan", data=stanDat, control = list(adapt_delta = 0.999))
  
  return(stanFit4=stanFit4)
}



# Function to get mode

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


# Choose 5k sites

set.seed(555)
indice <- sample(1:dim(FullAnnotation)[1],5000,replace=FALSE)
indice <- sort(indice)
indice_6 <- indice[2501:3000]
nsites <- length(indice_6)



mu_S <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                   p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                   p97.5=numeric(nsites))

sigmaE_S <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                       p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                       p97.5=numeric(nsites))


sigmaP_S <- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                       p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                       p97.5=numeric(nsites))


sigmaT_S<- data.frame(mean=numeric(nsites),mode=numeric(nsites),p2.5=numeric(nsites),
                      p25=numeric(nsites),p50=numeric(nsites),p75=numeric(nsites),
                      p97.5=numeric(nsites))



count <- 1

for (i in indice_6) {
  
  data <- site(i)
  stan <- stanfit4(data)
  posterior <- as.array(stan)
  
  
  mu_S[count,1] <- summary(stan)$summary[49,1]
  mu_S[count,2] <- getmode(posterior[1:1000,1:4,49])
  mu_S[count,3:7] <- summary(stan)$summary[49,4:8]
  
  
  sigmaE_S[count,1] <- summary(stan)$summary[50,1]
  sigmaE_S[count,2] <- getmode(posterior[1:1000,1:4,50])
  sigmaE_S[count,3:7] <- summary(stan)$summary[50,4:8]
  
  
  sigmaP_S[count,1] <- summary(stan)$summary[51,1]
  sigmaP_S[count,2] <- getmode(posterior[1:1000,1:4,51])
  sigmaP_S[count,3:7] <- summary(stan)$summary[51,4:8]
  
  
  sigmaT_S[count,1] <- summary(stan)$summary[52,1]
  sigmaT_S[count,2] <- getmode(posterior[1:1000,1:4,52])
  sigmaT_S[count,3:7] <- summary(stan)$summary[52,4:8]
  
  
  count <- count + 1
}

# save data
save(mu_S,sigmaP_S,sigmaT_S,sigmaE_S,file="STANS6.RData")
