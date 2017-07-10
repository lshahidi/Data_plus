# load StanCRun data
# produce full size data frames for mu,betaT,sigmas (P,PT,T,E)
# read in individual save files containing 1000 site runs
# store 1000 site chunks in corresponding positions in full data frame
# repeat for each completed run

library(ggplot2)
library(reshape2)

# Set working directory to normal data wd then folder StanCRuns_6K
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data/StanCRuns_6K")

betaT_Cfull <-
  data.frame(
    mean = numeric(866836),
    mode = numeric(866836),
    p2.5 = numeric(866836),
    p25 = numeric(866836),
    p50 = numeric(866836),
    p75 = numeric(866836),
    p97.5 = numeric(866836)
  )

mu_Cfull <- betaT_Cfull
sigmaE_Cfull <- betaT_Cfull
sigmaP_Cfull <- betaT_Cfull
sigmaPT_Cfull <- betaT_Cfull
sigmaT_Cfull <- betaT_Cfull

filesPresent <- c(1, 2, 3, 4, 6, 7)
fileNames <- dir(pattern = "StanCResults_*")
allInds <- numeric()
for (i in filesPresent) {
  load(fileNames[filesPresent == i])
  inds <- (1:1000) + (i - 1) * 1000
  allInds <- append(allInds, inds)
  betaT_Cfull[inds, ] <- betaT_C
  mu_Cfull[inds, ] <- mu_C
  sigmaE_Cfull[inds, ] <- sigmaE_C
  sigmaP_Cfull[inds, ] <- sigmaP_C
  sigmaPT_Cfull[inds, ] <- sigmaPT_C
  sigmaT_Cfull[inds, ] <- sigmaT_C
}


# plot variances
sigmaLMER2 <-
  data.frame(
    sigmaE = sigmaE_Cfull$p50[allInds],
    sigmaP = sigmaP_Cfull$p50[allInds],
    sigmaPT = sigmaPT_Cfull$p50[allInds],
    sigmaT = sigmaT_Cfull$p50[allInds]
  )
gg <- melt(sigmaLMER2)
ggplot(gg, aes(x = value, fill = variable)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(variable ~ .) +ggtitle("Medians of sigmas (6000 sites)") +xlim(0,2)

sigmaLMER2 <-
  data.frame(
    sigmaE = sigmaE_Cfull$mean[allInds],
    sigmaP = sigmaP_Cfull$mean[allInds],
    sigmaPT = sigmaPT_Cfull$mean[allInds],
    sigmaT = sigmaT_Cfull$mean[allInds]
  )
gg <- melt(sigmaLMER2)
ggplot(gg, aes(x = value, fill = variable)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(variable ~ .) +ggtitle("Means of sigmas (6000 sites)") +xlim(0,2)

sigmaLMER2 <-
  data.frame(
    sigmaE = sigmaE_Cfull$mode[allInds],
    sigmaP = sigmaP_Cfull$mode[allInds],
    sigmaPT = sigmaPT_Cfull$mode[allInds],
    sigmaT = sigmaT_Cfull$mode[allInds]
  )
gg <- melt(sigmaLMER2)
ggplot(gg, aes(x = value, fill = variable)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(variable ~ .) +ggtitle("Modes of sigmas (6000 sites)") +xlim(0,2)
