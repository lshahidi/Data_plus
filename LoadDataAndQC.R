#################################################################################
######## This program is used for reading in dataset and qualify control ########
######### Refer to the RMarkdown file (coming soon) for final results ###########
#################################################################################


library(minfi)
library(minfiData)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(gdata) #gdata package needs Perl language installed
library(RColorBrewer)

# here set the working directory that points to the data folder
# e.g. the folder with all datasets in it, should contain all the
# 1337-1387 folders
# please just comment out the directories already here, so each
# user can uncomment the setwd for them as we move code around

# Kevin's working directory
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# Yanlin's working directory
setwd("D:/DataPlus2017/Data")

############### Function used for read in data ####################

# base_dir = directory containing basenames file
# targets = name of the dataset after reading in basenames file
# pat_dir = directory containing patient information file
# pat_file = WHOLE name (including directory) of the dataset after reading in patient information file
# work_name = name of the RGChannelSet

read.fun <- function(base_dir,targets_name,work_name,pat_file,pat_name) {
  # read in patient data from the .xlsx file
  pat_name <- read.xls(pat_file,header=TRUE)
  
  # Extract targets
  targets_name <- data.frame(pat_name[,"Complete.Barcode"])
  colnames(targets_name) <- "Basement"
  
  # Read in all .idat files
  targets_name$Basement <- file.path(base_dir, targets_name$Basement)
  work_name <- read.metharray(targets_name$Basement, verbose = TRUE)
  
  # Push the RGChannelSet to the global environment
  return(work_name)
}


############################### Read in data  #################################

## 1337 ##

base_dir_1337 <- "1337_Shibata EPIC DNA methylation data package/IDAT FILES"
targets_name_1337 <- "targets_1337"
pat_file_1337 <- "1337_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1337 (Shibata-8).xls"
pat_name_1337 <- "pat_1337"
work_name_1337 <- "work_1337"

work_1337 <- read.fun(base_dir_1337,targets_name_1337,work_name_1337,pat_file_1337,pat_name_1337)



## 1345 ##

base_dir_1345 <- "1345_Shibata EPIC DNA methylation data package/IDAT FILES"
targets_name_1345 <- "targets_1345"
pat_file_1345 <- "1345_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1345 (Shibata-16).xlsx"
pat_name_1345 <- "pat_1345"
work_name_1345 <- "work_1345"

work_1345 <- read.fun(base_dir_1345,targets_name_1345,work_name_1345,pat_file_1345,pat_name_1345)



## 1350 ##

base_dir_1350 <- "1350_SHIBATA EPIC DNA METHYLATION DATA PACKAGE/IDAT files"
targets_name_1350 <- "targets_1350"
pat_file_1350 <- "1350_SHIBATA EPIC DNA METHYLATION DATA PACKAGE/SAMPLE-ARRAY MAPPING/1350 (Shibata-8).xlsx"
pat_name_1350 <- "pat_1350"
work_name_1350 <- "work_1350"

work_1350 <- read.fun(base_dir_1350,targets_name_1350,work_name_1350,pat_file_1350,pat_name_1350)
 

## 1357 ##

base_dir_1357 <- "1357_Shibata EPIC DNA methylation data package/IDAT FILES"
targets_name_1357 <- "targets_1357"
pat_file_1357 <- "1357_Shibata EPIC DNA Methylation Data Package/SAMPLE-ARRAY MAPPING/1357 (Shibata-8).xlsx"
pat_name_1357 <- "pat_1357"
work_name_1357 <- "work_1357"

work_1357 <- read.fun(base_dir_1357,targets_name_1357,work_name_1357,pat_file_1357,pat_name_1357)


## 1360 ##

base_dir_1360 <- "1360_Shibata EPIC Data Package/IDAT FILES"
targets_name_1360 <- "targets_1360"
pat_file_1360 <- "1360_Shibata EPIC Data Package/SAMPLE-ARRAY MAPPING/1360 (Shibata-8).xlsx"
pat_name_1360 <- "pat_1360"
work_name_1360 <- "work_1360"

work_1360 <- read.fun(base_dir_1360,targets_name_1360,work_name_1360,pat_file_1360,pat_name_1360)

## 1378 ##

base_dir_1378 <- "1378_Shibata EPIC Data Package/IDAT FILES"
targets_name_1385 <- "targets_1378"
pat_file_1378 <- "1378_Shibata EPIC Data Package/SAMPLE-ARRAY MAPPING/1378 (Shbata-8).xls"
pat_name_1378 <- "pat_1378"
work_name_1378 <- "work_1385"

work_1378 <- read.fun(base_dir_1378,targets_name_1378,work_name_1378,pat_file_1378,pat_name_1378)

## 1385 ##

base_dir_1385 <- "1385_Shibata EPIC Data Package/IDAT FILES"
targets_name_1385 <- "targets_1385"
pat_file_1385 <- "1385_Shibata EPIC Data Package/SAMPLE-ARRAY MAPPING/1385 (Shibata-8).xlsx"
pat_name_1385 <- "pat_1385"
work_name_1385 <- "work_1385"

work_1385 <- read.fun(base_dir_1385,targets_name_1385,work_name_1385,pat_file_1385,pat_name_1385)

 
## 1387 ##

base_dir_1387 <- "1387_Shibata EPIC DNA methylation data package/IDAT FILES"
targets_name_1387 <- "targets_1387"
pat_file_1387 <- "1387_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1387 (Shibata-8).xls"
pat_name_1387 <- "pat_1387"
work_name_1387 <- "work_1387"

work_1387 <- read.fun(base_dir_1387,targets_name_1387,work_name_1387,pat_file_1387,pat_name_1387)



################################## Bad Sample in 1337 ####################################

# Convert to MethylSet
mset_1337 <- preprocessRaw(work_1337)

# Check the methylation and unmethylation signals
head(getMeth(mset_1337))
head(getUnmeth((mset_1337)))

# Convert to RatioSet containing beta values and Mvalues(log(M/U))
rset_1337 <- ratioConvert(mset_1337,what = "both")

# Get beta
beta_1337 <- getBeta(rset_1337)

# Check beta values
head(beta_1337)

#200360140022_R07C01 (JB) has bad quality
dim(beta_1337)[1] == sum(is.na(beta_1337[,7]))

### read in new data without JB ###
base_dir_1337_new <- "1337_Shibata EPIC DNA methylation data package/IDAT FILES"

targets_1337_new <- read.csv("1337_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1337_Shibata-targets.csv", as.is = TRUE)

targets_1337_new$Basename <- file.path(base_dir_1337_new, targets_1337_new$Basename)
work_1337_new <- read.metharray(targets_1337_new$Basename, verbose = TRUE)



############################## Preprocessing/Normalization ###############################

### 1337 ###

## noob ##
noob_1337 <- preprocessNoob(work_1337_new)

rset_1337_noob <- ratioConvert(noob_1337,what="both")
beta_1337_noob <- getBeta(rset_1337_noob)
head(beta_1337_noob)
# Compare with raw beta values
head(beta_1337)


## Genome Studio ##
lumi_1337 <- preprocessIllumina(work_1337_new,bg.correct = TRUE,
                                normalize = "controls",reference = 2)
rset_1337_lumi <- ratioConvert(lumi_1337,what="both")
beta_1337_lumi <- getBeta(rset_1337_lumi)
head(beta_1337_lumi)
# Compare with noob and raw
head(beta_1337)
head(beta_1337_noob)


### 1345 ###

## noob ##
noob_1345 <- preprocessNoob(work_1345)

## Genome Studio ##
lumi_1345 <- preprocessIllumina(work_1345,bg.correct = TRUE,
                                normalize = "controls",reference = 2)


### 1350 ###

## noob ##
noob_1350 <- preprocessNoob(work_1350)

## Genome Studio ##
lumi_1350 <- preprocessIllumina(work_1350,bg.correct = TRUE,
                                normalize = "controls",reference = 2)


### 1357 ###

## noob ##
noob_1357 <- preprocessNoob(work_1357)

## Genome Studio ##
lumi_1357 <- preprocessIllumina(work_1357,bg.correct = TRUE,
                                normalize = "controls",reference = 2)


### 1360 ###

## noob ##
noob_1360 <- preprocessNoob(work_1360)

## Genome Studio ##
lumi_1360 <- preprocessIllumina(work_1360,bg.correct = TRUE,
                                normalize = "controls",reference = 2)


### 1378 ###

## noob ##
noob_1378 <- preprocessNoob(work_1378)

## Genome Studio ##
lumi_1378 <- preprocessIllumina(work_1378,bg.correct = TRUE,
                                normalize = "controls",reference = 2)


### 1385 ###

## noob ##
noob_1385 <- preprocessNoob(work_1385)

## Genome Studio ##
lumi_1385 <- preprocessIllumina(work_1385,bg.correct = TRUE,
                                normalize = "controls",reference = 2)


### 1387 ###

## noob ##
noob_1387 <- preprocessNoob(work_1387)

## Genome Studio ##
lumi_1387 <- preprocessIllumina(work_1387,bg.correct = TRUE,
                                normalize = "controls",reference = 2)



#################################### Quality Control #####################################

# Beta density plot
pat <- read.xls("1337_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1337 (Shibata-8).xls")

par(mfrow=c(1,2), mar=c(4,2,4,1))
densityPlot(noob_1337,sampGroups=pat$Sample_No,main = "Corrected Beta Values Using 'noob'",
            legend = FALSE)
legend("topright",legend = levels(pat$Sample_No),text.col=brewer.pal(8,"Dark2"),cex=0.7)
densityPlot(lumi_1337,sampGroups=pat$Sample_No,main = "Corrected Beta Values Using 'Illumina'",
            legend = FALSE)
legend("topright",legend = levels(pat$Sample_No),text.col=brewer.pal(8,"Dark2"),cex=0.7)

qc_1337_noob <- getQC(noob_1337)
plotQC(qc_1337_noob)
qc_1337_lumi <- getQC(lumi_1337)
plotQC(qc_1337_lumi)


# Plot sex
gmset_1337 <- mapToGenome(noob_1337)
sex_1337 <- getSex(gmset_1337)
plotSex(sex_1337)

# MDS plot
mdsPlot(noob_1337)

qcReport(work_1337_new,sampNames = pat$Sample_No[-7],
         sampGroups = pat$Plate[-7],pdf = "qcReport.pdf")

# this line works to show density curves :)
densityPlot(work_1337_new, sampGroups = pat$Sample_No, main = "Beta", xlab = "Beta")

# insert more QC


## Pre-processing

# built-in preprocessing, we should do this ourselves
mset <- preprocessIllumina(RGset)
mset <- mapToGenome(mset)

dim(getBeta(mset, type = "Illumina"))  ##the argument type='Illumina' gives us default procedure
head(granges(mset))

sex <- getSex(mset)
plotSex(sex)

plot(as.matrix(getQC(mset)))
