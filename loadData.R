#################################################################################
######## This program is used for reading in dataset and qualify control ########
######### Refer to the RMarkdown file (coming soon) for final results ###########
#################################################################################


library(minfi)
library(minfiData)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(gdata)

# here set the working directory that points to the data folder
# e.g. the folder with all datasets in it, should contain all the
# 1337-1387 folders
# please just comment out the directories already here, so each
# user can uncomment the setwd for them as we move code around
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")
# setwd("put yours here")

############### Function used for read in data ####################

# base_dir = directory containing basenames file
# targets_name = name of the dataset after reading in basenames file
# pat_dir = directory containing patient information file
# pat_file = WHOLE name (including directory) of the dataset after reading in patient information file
# work_name = name of the RGChannelSet

read.fun <- function(base_dir,targets_name,work_name,pat_file,pat_name) {
  # read in patient data from the .xls file
  pat_name <- read.xls(pat_file, sheetIndex=1,header=TRUE)
  
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

# Kevin's version
base_dir <- "/Users/kevinmurgas/Documents/Data+ project/EPIC data/1337_Shibata EPIC DNA methylation data package/IDAT FILES"
targets_name <- "targets_1337"
pat_file <- "/Users/kevinmurgas/Documents/Data+ project/EPIC data/1337_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1337 (Shibata-8).xls"
pat_name <- "pat_1337"
work_name <- "work_1337"

#work_name you chose# <- read.fun(base_dir,targets_name,work_name,pat_file,pat_name)

# Yanlin's version

base_dir_1337 <- "D:/DataPlus2017/Data/1337_Shibata EPIC DNA methylation data package/1337_Shibata EPIC DNA methylation data package/IDAT FILES"
targets_name_1337 <- "targets_1337"
pat_file_1337 <- "D:/DataPlus2017/Data/1337_Shibata EPIC DNA methylation data package/1337_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1337 (Shibata-8).xls"
pat_name_1337 <- "pat_1337"
work_name_1337 <- "work_1337"

work_1337 <- read.fun(base_dir_1337,targets_name_1337,work_name_1337,pat_file_1337,pat_name_1337)



## 1345 ##

# Kevin's version

# Yanlin's version
base_dir_1345 <- "D:/DataPlus2017/Data/1345_Shibata EPIC DNA methylation data package/1345_Shibata EPIC DNA methylation data package/IDAT FILES"
targets_name_1345 <- "targets_1345"
pat_file_1345 <- "D:/DataPlus2017/Data/1345_Shibata EPIC DNA methylation data package/1345_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1345 (Shibata-16).xlsx"
pat_name_1345 <- "pat_1345"
work_name_1345 <- "work_1345"

work_1345 <- read.fun(base_dir_1345,targets_name_1345,work_name_1345,pat_file_1345,pat_name_1345)



## 1350 ##

# Yanlin' version
base_dir_1350 <- "D:/DataPlus2017/Data/1350_SHIBATA EPIC DNA METHYLATION DATA PACKAGE/1350_SHIBATA EPIC DNA METHYLATION DATA PACKAGE/IDAT files"
targets_name_1350 <- "targets_1350"
pat_file_1350 <- "D:/DataPlus2017/Data/1350_SHIBATA EPIC DNA METHYLATION DATA PACKAGE/1350_SHIBATA EPIC DNA METHYLATION DATA PACKAGE/SAMPLE-ARRAY MAPPING/1350 (Shibata-8).xlsx"
pat_name_1350 <- "pat_1350"
work_name_1350 <- "work_1350"

work_1350 <- read.fun(base_dir_1350,targets_name_1350,work_name_1350,pat_file_1350,pat_name_1350)
 

## 1357 ##

# Yanlin's version
base_dir <- "D:/DataPlus2017/Data/1357_Shibata EPIC DNA methylation data package/1357_Shibata EPIC DNA methylation data package/IDAT FILES"
targets_name <- "targets_1357"
pat_file <- "D:/DataPlus2017/Data/1357_Shibata EPIC DNA methylation data package/1357_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1357 (Shibata-8).xlsx"
pat_name <- "pat_1357"
work_name <- "work_1357"

work_1357 <- read.fun(base_dir,targets_name,work_name,pat_file,pat_name)


## 1360 ##

# Yanlin's version
base_dir <- "D:/DataPlus2017/Data/1360_Shibata EPIC DNA methylation data package/1360_Shibata EPIC DNA methylation data package/IDAT FILES"
targets_name <- "targets_1345"
pat_file <- "D:/DataPlus2017/Data/1345_Shibata EPIC DNA methylation data package/1345_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1345 (Shibata-16).xlsx"
pat_name <- "pat_1345"
work_name <- "work_1345"

work_1345 <- read.fun(base_dir,targets_name,work_name,pat_file,pat_name)

## 1378 ##

## 1385 ##
 
## 1387 ##



#################################### Old version of reading in data #############################################
############################### you can delete if the function above works ########################################

# here set a working directory in the folder with EPIC data
# specifically the directory containing folders for 1337, 1345, etc
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# start by choosing 1337 as base directory
baseDir <- "1337_Shibata EPIC DNA methylation data package/IDAT FILES"

# read in BaseNames from the .csv file
targets <- read.csv("1337_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1337_Shibata-targets.csv", as.is = TRUE)

# converts target files to RGset
targets$Basename <- file.path(baseDir, targets$Basename)
RGset <- read.metharray(targets$Basename, verbose = TRUE)
# annotation(RGsetEx) # check annotation packages

# read in patient data from the .xls file
pd <- read.xls("1337_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1337 (Shibata-8).xls", sheet = 1, header = TRUE)
####################################################################################################################






#################################### Quality Control #####################################

qcReport(work_1337)


# vvv this line does not work ¯\_(ツ)_/¯
# qcReport(RGset, sampNames = pd$Sample_No, sampGroups = pd$Plate, pdf = "qcReport.pdf")

# this line works to show density curves :)
densityPlot(RGset, sampGroups = pd$Sample_No, main = "Beta", xlab = "Beta")

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
