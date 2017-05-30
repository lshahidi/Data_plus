#################################################################################
######## This program is used for reading in dataset and qualify control ########
######### Refer to the RMarkdown file (coming soon) for final results ###########
#################################################################################


library(minfi)
library(minfiData)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(gdata)


############### Function used for read in data ####################

# base_dir = directory containing basenames file
# targets_name = name of the dataset after reading in basenames file
# pat_dir = directory containing patient information file
# pat_file = WHOLE name (including directory) of the dataset after reading in patient information file
# work_name = name of the RGChannelSet

read.fun <- function(base_dir,targets_name,work_name,pat_file,pat_name) {
  # read in patient data from the .xls file
  pat_name <- read.xlsx(pat_file, sheetIndex=1,header=TRUE)
  
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
targets_name <- #choose what you like#
pat_file <- "/Users/kevinmurgas/Documents/Data+ project/EPIC data/1337_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1337 (Shibata-8).xls"
pat_name <- #choose what you like#
work_name <- #choose what you like#

#work_name you chose# <- read.fun(base_dir,targets_name,work_name,pat_file,pat_name)

# Yanlin's version

base_dir <- "D:/DataPlus2017/Data/1337_Shibata EPIC DNA methylation data package/1337_Shibata EPIC DNA methylation data package/IDAT FILES"
targets_name <- "targets_1337"
pat_file <- "D:/DataPlus2017/Data/1337_Shibata EPIC DNA methylation data package/1337_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1337 (Shibata-8).xls"
pat_name <- "pat_1337"
work_name <- "work_1337"

work_1337 <- read.fun(base_dir,targets_name,work_name,pat_file,pat_name)



## 1345 ##

# Kevin's version

# Yanlin's version
base_dir <- "D:/DataPlus2017/Data/1345_Shibata EPIC DNA methylation data package/1345_Shibata EPIC DNA methylation data package/IDAT FILES"
targets_name <- "targets_1345"
pat_file <- "D:/DataPlus2017/Data/1345_Shibata EPIC DNA methylation data package/1345_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1345 (Shibata-16).xlsx"
pat_name <- "pat_1345"
work_name <- "work_1345"

work_1345 <- read.fun(base_dir,targets_name,work_name,pat_file,pat_name)



## 1350 ##
 
## 1357 ##

## 1360 ##

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
