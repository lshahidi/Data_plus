# script for extracting variance from Stan model 
# part 2: on second half


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
  stanFit2 <- stan(file="model2.stan", data=stanDat)

  return(stanFit2=stanFit2)
}


nsites <- dim(FullAnnotation)[1]/2

Data <- list(NULL)

for (i in (nsites+1):(2*nsites)) {
  Data[[count]] <- site(i)
}

stan <- lapply(Data, stanfit)

sigmaSTAN <- as.data.frame(matrix(NA,nsites,3))


for (i in 1:nsites) {
  sigmaSTAN[i,] <- summary(stan[[i]])$summary[35:37,1]
}

colnames(sigmaSTAN) <- c("sigmaE","sigmaP","sigmaT")
rownames(sigmaSTAN) <- paste("site",(nsites+1):(2*nsites))

sigmaSTAN <- sigmaSTAN[,c(3,2,1)]

# save data
save(sigmaSTAN,file="mySTANVars2.RData")
