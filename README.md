# Quantifying Tumor Evolution with Hierarchical Statistical Modeling

This repository contains code used to work with 850K EPIC methylation array data in order to develop a model for variation of methylation among multiple patients and normal/tumor tissues.

README updated: 08/08/2017 KM

## Files Present

Scripts: .R files used to run code

* compareLMER_STANmodels.R - compare between lmer and Stan results
* compareLMERmodels.R - used to compare between two lmer models (correlating and non-correlating)
* compareSTANmodels.R - outdated, used to compare between two different Stan models (model.stan and model2.stan)
* DataModel.R - testing script, sort of bucket of code, contains lmer model to produce small quantities of data, for plotting results, defining gene information, scoring genes, comparing results
* individualGenes.R - script for examining individual genes, including scoring and correlation plots
* lmerVarsScript1.R - outdated, script to run lmer model on 1000 sites
* LoadDataAndQC.R - load raw EPIC array data from 8 provided data sets, preprocess with noob, perform QC, and annotate. Saves as "myFA.Rdata"
* loadSCRuns.R - loads and combines results from parallel runs
* NewStanCParallel.R - similar to StanCParallel.R, runs Stan model using model3.stan code, including PTprob statistic (sum of posterior probability for: log(sigmaP/sigmaT) > 0 )
* Score.R - testing code used to develop gene scoring methods
* ScoreGenesParallel.R - script to score genes in parallel
* StanCParallel.R - runs Stan model using model3.stan code, to extract at each site of 10000 sites chunk the following variables: mu, betaT, sigmaE, sigmaP, sigmaPT, sigmaT, then save  in a numbered result file
* StanCParallel2.R - testing code for parallelization, only runs 3 sites, no save
* StanCVarsScript1.R - outdated, used to run model3.stan on 5000 sites

Markdown files: .Rmd files used to produce markdown documents with R code

* lmer Model.Rmd - runs three unique lmer models on all data and presents results
* Model.Rmd - runs one lmer model and presents results
* QC_and_Annotation.Rmd - loads data and presents quality control results
* stan_Model.Rmd - runs Stan model using model.stan (out of date)

Stan Models: .stan files used to define Stan model equations

* model.stan - original model, with parameters mu, betaT, b and bT
* model2.stan - same as model.stan but now including priors for mu and betaT
* model3.stan - complex model, includes normal/tumor indicator and three levels of variance (P, PT, T) with parameters mu, betaT, b, c, d, 
* model4.stan - simplified model, without normal/tumor indicator, not in use
* model_rep.stan - model3.stan reparameterized, never finalized

Misc files:

* Data+ 2017.Rproj - contains data for RStudio project
* README.md - c'est moi

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
  * Score.R - scores genes using 5 unique methods
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
* gdata
* ggplot2
* gtools
* IlluminaHumanMethylationEPICmanifest
* IlluminaHumanMethylationEPICanno.ilm10b2.hg19
* minfi
* minfiData
* reshape2
* RColorBrewer
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
