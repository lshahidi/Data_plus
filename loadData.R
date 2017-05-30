library(minfi)
library(minfiData)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# here set a working directory in the folder with EPIC data
# specifically the directory containing folders for 1337, 1345, etc
setwd("/Users/kevinmurgas/Documents/Data+ project/EPIC data")

# start by choosing 1337
baseDir <- "1337_Shibata EPIC DNA methylation data package/IDAT FILES"

targets <- read.csv("1337_Shibata EPIC DNA methylation data package/SAMPLE-ARRAY MAPPING/1337_Shibata-targets.csv", as.is = TRUE)

sub(baseDir, "", targets$Basename)

targets$Basename <- file.path(path, targets$Basename)
RGset <- read.metharray(targets$Basename, verbose = TRUE)

# insert quality control

# built-in preprocessing, we should do this ourselves
mset <- preprocessIllumina(RGset)
mset <- mapToGenome(mset)

dim(getBeta(mset, type = "Illumina"))  ##the argument type='Illumina' gives us default procedure


