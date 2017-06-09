### plot beta values between specific samples

### INITIALIZE

library(ggplot2)

# here set the working directory that points to the data folder
# e.g. the folder with annotated data saved as "myFA.Rdata"
# please just comment out the directories already here, so each
# user can uncomment the setwd for them as we move code around

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
setwd("D:/DataPlus2017/Data")


##### LOAD DATA

# load fully annotated data (saved from LoadDataAndQC.R)
load("myFA.Rdata")


##### COMPARE BETA

### tumor vs normal for patients H E K* W C J (have both)
# average tumor side A and B for tumor data

# H
betaHNor <- FullAnnotation$HN
betaHTum <- (FullAnnotation$HA + FullAnnotation$HB)/2
Hdata <- data.frame(x=betaHNor, y=betaHTum)
ggplot(Hdata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta H Normal", y="Beta H Tumor")
cov(betaHNor, betaHTum)
c(mean(betaHNor), mean(betaHTum), mean(betaHNor-betaHTum))

# E
betaENor <- FullAnnotation$EN
betaETum <- (FullAnnotation$EA + FullAnnotation$EB)/2
Edata <- data.frame(x=betaENor, y=betaETum)
ggplot(Edata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta E Normal", y="Beta E Tumor")
cov(betaENor, betaETum)
c(mean(betaENor), mean(betaETum), mean(betaENor-betaETum))

# K*
betaKnNor <- FullAnnotation$`KN (new)`
betaKnTum <- (FullAnnotation$`KA (new)` + FullAnnotation$`KB (new)`)/2
Kndata <- data.frame(x=betaKnNor, y=betaKnTum)
ggplot(Kndata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta K* Normal", y="Beta K* Tumor")
cov(betaKnNor, betaKnTum)
c(mean(betaKnNor), mean(betaKnTum), mean(betaKnNor-betaKnTum))

# W
betaWNor <- FullAnnotation$WN
betaWTum <- (FullAnnotation$WA + FullAnnotation$WB)/2
Wdata <- data.frame(x=betaWNor, y=betaWTum)
ggplot(Wdata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta W Normal", y="Beta W Tumor")
cov(betaWNor, betaWTum)
c(mean(betaWNor), mean(betaWTum), mean(betaWNor-betaWTum))

# C
betaCNor <- FullAnnotation$CN
betaCTum <- (FullAnnotation$CA + FullAnnotation$CB)/2
Cdata <- data.frame(x=betaCNor, y=betaCTum)
ggplot(Cdata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta C Normal", y="Beta C Tumor")
cov(betaCNor, betaCTum)
c(mean(betaCNor), mean(betaCTum), mean(betaCNor-betaCTum))

# J
betaJNor <- FullAnnotation$JN
betaJTum <- (FullAnnotation[,39] + FullAnnotation$JB)/2
Jdata <- data.frame(x=betaJNor, y=betaJTum)
ggplot(Jdata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta J Normal", y="Beta J Tumor")
cov(betaJNor, betaJTum)
c(mean(betaJNor), mean(betaJTum), mean(betaJNor-betaJTum))

### tumor side A vs side B

# K
betaKA <- FullAnnotation$KA
betaKB <- FullAnnotation$KB
Ktumor <- data.frame(x=betaKA, y=betaKB)
ggplot(Ktumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta K Side A", y="Beta K Side B")
cov(betaKA, betaKB)

# H
betaHA <- FullAnnotation$HA
betaHB <- FullAnnotation$HB
Htumor <- data.frame(x=betaHA, y=betaHB)
ggplot(Htumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta H Side A", y="Beta H Side B")
cov(betaHA, betaHB)

# E
betaEA <- FullAnnotation$EA
betaEB <- FullAnnotation$EB
Etumor <- data.frame(x=betaEA, y=betaEB)
ggplot(Etumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta E Side A", y="Beta E Side B")
cov(betaEA, betaEB)

# K*
betaKnA <- FullAnnotation$`KA (new)`
betaKnB <- FullAnnotation$`KB (new)`
Kntumor <- data.frame(x=betaKnA, y=betaKnB)
ggplot(Kntumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta K* Side A", y="Beta K* Side B")
cov(betaKnA, betaKnB)

# X
betaXA <- FullAnnotation$XA
betaXB <- FullAnnotation$XB
Xtumor <- data.frame(x=betaXA, y=betaXB)
ggplot(Xtumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta X Side A", y="Beta X Side B")
cov(betaXA, betaXB)

# W
betaWA <- FullAnnotation$WA
betaWB <- FullAnnotation$WB
Wtumor <- data.frame(x=betaWA, y=betaWB)
ggplot(Wtumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta W Side A", y="Beta W Side B")
cov(betaWA, betaWB)

# T
betaTA <- FullAnnotation$TA
betaTB <- FullAnnotation$TB
Ttumor <- data.frame(x=betaTA, y=betaTB)
ggplot(Ttumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta T Side A", y="Beta T Side B")
cov(betaTA, betaTB)

# S
betaSA <- FullAnnotation$SA
betaSB <- FullAnnotation$SB
Stumor <- data.frame(x=betaSA, y=betaSB)
ggplot(Stumor,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta S Side A", y="Beta S Side B")
cov(betaSA, betaSB)
mean(betaSA)-mean(betaSB)
