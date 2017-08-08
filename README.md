# Quantifying Tumor Evolution with Hierarchical Statistical Modeling

This repository contains code used to work with 850K EPIC methylation array data in order to develop a model for variation of methylation among multiple patients and normal/tumor tissues.

## Files Present

Scripts: .R files used to run code

* compareLMER_STANmodels.R - 
* compareLMERmodels.R -
* compareSTANmodels.R -
* DataModel.R - 
* individualGenes.R - 
* lmerVarsScript1.R - 
* LoadDataAndQC.R

Markdown files: .Rmd files used to produce markdown documents with R code

* Hat tip to anyone who's code was used
* Inspiration
* etc

Scripts: .R files used to run code

* Hat tip to anyone who's code was used
* Inspiration
* etc

Misc files:

* Hat tip to anyone who's code was used
* Inspiration
* etc

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

One must obtain the EPIC data sampled by Dr. Daryl Shibata in order to produce the data sets that this code was developed to work with.

Having the correct data, one can use the following workflow:

1. *Load Data, Preprocess and QC* - Run **LoadDataAndQC.R** to process raw data to normalized beta values and perform quality control.
  * to begin simply add the working directory that points to the folder containing EPIC data (e.g. add in 'setwd("YOURWORKINGDIRECTORY")' )
  * please leave previous users' working directory lines for convenience
  * this script also adds data annotation read from a .csv file
  
2. *Model Data with Stan* - With the beta values preprocessed, the next step is to perform hierarchical statistical modeling. We have the choice between several lmer and RStan models. For our project purposes, we decided on RStan using model3.stan. We also *highly* recommend using parallel processing methods
  * Stan scripts:
    * **StanCParallel.R** - this is the parallelized script designed for use on the DukeComputerCluster, using RStan, model3.stan (C='Complex model') on 10000 site chunks. Takes ~2 days to run
    * **NewStanCParallel.R** - same as StanCParallel.R, now including PTprob statistic
    * StanCParallel2.R - this script used only for testing StanCParallel.R code on local machine (modified for input, reduced sites run, and no saving)
    * StanCVarsScript1.R - outdated, preceded StanCParallel.R before parallelization and optimization
  * lmer Scripts:
    * lmerVarsScript1.R - runs lmer model (contained in code) in 1000 site chunks
    
2.5. *Load Parallel Runs* - when using parallel methods on a computer cluster, the data is exported in chunks and thus must be loaded and combined with **loadSCRuns.R**
    
3. *Analyze Results* - we can run a number of scripts to look at lmer and RStan results
  * compareSTANmodels.R - outdated, used to compare between two different Stan models (model.stan and model2.stan)
  * compareLMER_STANmodels.R - compare between lmer and Stan results
  * compareLMERmodels.R - used to compare between two lmer models (correlating and non-correlating)
  
4. *Score Genes* - run these scripts to test gene scoring methods, given model results (StanC)
  * Score.R - scores genes using [] methods
  * ScoreGenesParallel.R - same script just modified to run in parallel on DukeComputerCluster, cuts time down to <10 minutes

### Prerequisites

The following packages are used throughout the scripts, thus must be installed manually via:

```
install.packages("PACKAGENAME")
```

*Packages*:
* arm
* bayesplot
* coda
* doParallel
* lme4
* ggplot2
* gtools
* reshape2
* rstan

## Authors

* **Yanlin Ma** - *Student Intern* - equal contributions to the project and coding
* **Kevin Murgas** - *Student Intern* - equal contributions to the project and coding

Additional Contributions:
* **Marc D. Ryser, Ph.D.** - *Project Lead*
* **Daryl Shibata, M.D.** - *Project Lead*
* **Lidea Shahidi** - *Project Manager*
* **Andrew Allen, Ph.D.** - *Collaboration*

## Acknowledgments

* we thank people for the help :)
