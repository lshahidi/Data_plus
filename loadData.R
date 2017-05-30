library(minfi)
library(minfiData)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(gdata)


## Read in data

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


## Quality Control

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
