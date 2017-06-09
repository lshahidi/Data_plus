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


### LOAD DATA

# load fully annotated data (saved from LoadDataAndQC.R)
load("myFA.Rdata")

# examine tumor vs normal for patients H E K* W C (have both)
# average tumor side A and B

# H
betaHNor <- FullAnnotation$HN
betaHTum <- (FullAnnotation$HA + FullAnnotation$HB)/2
Hdata <- data.frame(x=betaHNor, y=betaHTum)
ggplot(Hdata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta H Normal", y="Beta H Tumor")
cov(betaHNor, betaHTum)

# E
betaENor <- FullAnnotation$EN
betaETum <- (FullAnnotation$EA + FullAnnotation$EB)/2
Edata <- data.frame(x=betaENor, y=betaETum)
ggplot(Edata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta E Normal", y="Beta E Tumor")
cov(betaENor, betaETum)

# K*
betaKnNor <- FullAnnotation$`KN (new)`
betaKnTum <- (FullAnnotation$`KA (new)` + FullAnnotation$`KB (new)`)/2
Kndata <- data.frame(x=betaKnNor, y=betaKnTum)
ggplot(Kndata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta K* Normal", y="Beta K* Tumor")
cov(betaKnNor, betaKnTum)

# W
betaWNor <- FullAnnotation$WN
betaWTum <- (FullAnnotation$WA + FullAnnotation$WB)/2
Wdata <- data.frame(x=betaWNor, y=betaWTum)
ggplot(Wdata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta W Normal", y="Beta W Tumor")
cov(betaWNor, betaWTum)

# C
betaCNor <- FullAnnotation$CN
betaCTum <- (FullAnnotation$CA + FullAnnotation$CB)/2
Cdata <- data.frame(x=betaCNor, y=betaCTum)
ggplot(Cdata,aes(x=x,y=y)) + geom_point(alpha = 0.01) + labs(x="Beta C Normal", y="Beta C Tumor")
cov(betaCNor, betaCTum)
