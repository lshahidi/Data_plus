# script for extracting variance from Stan complex model (model3)
# input taskID corresponding to 10000 site chunks starting at site (N-1)*10000+1
# Now incorporating parallel processing with parallel package

library(doParallel)
library(rstan)
library(gtools)

rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())
print(paste("Cores: ", parallel::detectCores()))

# load data
load("myFA.Rdata")

args <- commandArgs(trailingOnly = TRUE)
print(args)
N <- as.numeric(args[1])
print(paste("Task #: ", N))
out.file <- args[2]

# Choose 10000 sites by N, and select chunk within FullAnnotation
siteInds <- (1:10000) + (N - 1) * 10000
if (N > 86) {
  siteInds <- siteInds[1:6836]
}
nsites <- length(siteInds)

print(paste("Making chunk of nsites =", nsites))
chunk <- FullAnnotation[siteInds,]


### FXNS ###
# Function used to read in data from each site
# changed to now use a "chunk" instead of FullAnnotation, for speed and memory
# requires index 1-1000 instead of siteIndex
site <- function (ind) {
  # extract CpG site xx to start
  temp <- chunk[ind,]
  
  # here we use all patient samples, excluding glands
  indices <- c(9:13, 15:39, 46, 47, 57, 58, 72:75)
  patientLabel <- substr(colnames(temp[indices]), 1, 1)
  patientLabel[10:12] <- "K*"
  sideLabel <- substr(colnames(temp[indices]), 2, 2)
  tissueLabel <- sideLabel
  tissueLabel[tissueLabel %in% c("A", "B")] <-
    "T"   #replace A and B with T
  tumorIndicator <- 1 * (tissueLabel == "T")
  
  work <-
    data.frame(
      beta = logit(t(temp[indices])),
      patient = patientLabel,
      tissue = tissueLabel,
      side = sideLabel,
      tInd = tumorIndicator
    )
  
  return(work)
}

# Function to build stan model for each site
stanfit3 <- function (dataset) {
  stanDat <- list(
    pID = as.integer(factor(dataset$patient)),
    tInd = dataset$tInd,
    N = nrow(dataset),
    P = nlevels(dataset$patient),
    y = dataset[, 1]
  )
  
  # Using Model 3: add intra-tumoral variances
  stanFit3 <-
    stan(
      file = "model3.stan",
      data = stanDat,
      control = list(adapt_delta = 0.999),
      refresh = 0
    )
  
  return(stanFit3 = stanFit3)
}

# Function for combining parallel runs of each fixef/sigma
comb <- function(x, ...) {
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}


### CODE ###

print("Generating data.frames")

betaT_C <-
  data.frame(
    mean = numeric(nsites),
    p2.5 = numeric(nsites),
    p25 = numeric(nsites),
    p50 = numeric(nsites),
    p75 = numeric(nsites),
    p97.5 = numeric(nsites)
  )

mu_C <-
  data.frame(
    mean = numeric(nsites),
    p2.5 = numeric(nsites),
    p25 = numeric(nsites),
    p50 = numeric(nsites),
    p75 = numeric(nsites),
    p97.5 = numeric(nsites)
  )

sigmaE_C <-
  data.frame(
    mean = numeric(nsites),
    p2.5 = numeric(nsites),
    p25 = numeric(nsites),
    p50 = numeric(nsites),
    p75 = numeric(nsites),
    p97.5 = numeric(nsites)
  )

sigmaP_C <-
  data.frame(
    mean = numeric(nsites),
    p2.5 = numeric(nsites),
    p25 = numeric(nsites),
    p50 = numeric(nsites),
    p75 = numeric(nsites),
    p97.5 = numeric(nsites)
  )

sigmaT_C <-
  data.frame(
    mean = numeric(nsites),
    p2.5 = numeric(nsites),
    p25 = numeric(nsites),
    p50 = numeric(nsites),
    p75 = numeric(nsites),
    p97.5 = numeric(nsites)
  )

sigmaPT_C <-
  data.frame(
    mean = numeric(nsites),
    p2.5 = numeric(nsites),
    p25 = numeric(nsites),
    p50 = numeric(nsites),
    p75 = numeric(nsites),
    p97.5 = numeric(nsites)
  )

# parallel via doParallel
mInd <- c(1, 4:8)
registerDoParallel(detectCores())
print(paste("Cores registered:",getDoParWorkers()))
print("Starting foreach loop")
ptm <- proc.time()
parData <- foreach(i = (1:nsites), .combine = 'comb', .multicombine = TRUE) %dopar% {
  print(paste("site:", i))
  
  data <- site(i)
  stanFit <- stanfit3(data)
  fitSumm <- summary(stanFit)$summary[71:76, mInd]
  split(fitSumm, row(fitSumm))
}
proc.time() - ptm

betaT_C[,] <- parData$'1'
mu_C[,] <- parData$'2'
sigmaE_C[,] <- parData$'3'
sigmaP_C[,] <- parData$'4'
sigmaPT_C[,] <- parData$'5'
sigmaT_C[,] <- parData$'6'

print(paste("Completed run, now saving"))
# save data
save(mu_C, betaT_C, sigmaP_C, sigmaT_C,
     sigmaPT_C, sigmaE_C, file = out.file)