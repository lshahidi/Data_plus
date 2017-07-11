# script for extracting variance from Stan complex model (model3)
# input taskID corresponding to 1000 site chunks up to  N*1000  sites
# Now incorporating parallel processing with parallel package

#library(doParallel)
library(rstan)
library(gtools)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
print(paste("Cores: ",parallel::detectCores()))

# load data
load("myFA.Rdata")

args <- commandArgs(trailingOnly = TRUE)
print(args)
N <- as.numeric(args[1])
print(paste("Task #: ", N))
out.file <- args[2]

### FXNS ###
# Function used to read in data from each site
site <- function (site_no) {
  # extract CpG site xx to start
  temp <- FullAnnotation[site_no, ]
  
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

### CODE ###
# Choose 1000  sites by N

siteInds <- (1:1000) + (N - 1) * 1000
if (N > 866) {
  siteInds <- siteInds[1:836]
}
nsites <- length(siteInds)

print(paste("nsites: ", nsites))
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

#registerDoParallel(detectCores())
#getDoParWorkers()
#parData <- foreach(i = (1:nsites), .packages = c("rstan")) %dopar% {
mInd <- c(1,4:8)
ptm <- proc.time()
for (i in (1:nsites)) {
  print(paste("site: ",i))
  
  data <- site(siteInds[i])
  stanFit <- stanfit3(data)
  fitSumm <- summary(stanFit)
  
  betaT_C[i,] <- fitSumm$summary[71, mInd]
  mu_C[i,] <- fitSumm$summary[72, mInd]
  sigmaE_C[i,] <- fitSumm$summary[73, mInd]
  sigmaP_C[i,] <- fitSumm$summary[74, mInd]
  sigmaPT_C[i,] <- fitSumm$summary[75, mInd]
  sigmaT_C[i,] <- fitSumm$summary[76, mInd]
  
}
proc.time() - ptm

print(paste("Completed run, now saving"))
# save data
save(mu_C, betaT_C, sigmaP_C, sigmaT_C,
     sigmaPT_C, sigmaE_C, file = out.file)