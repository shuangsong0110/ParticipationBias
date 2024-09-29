# ParticipationBias

A pipeline for adjusting participation bias in the estimation of heritability and genetic correlation.


## Table of contents
* [Prerequisites](#white_check_mark-prerequisites)
* [Installation](#hammer_and_wrench-installation)
* [Prepare GWAS summary statistics](#scroll-prepare-gwas-summary-statistics)
* [Heritability adjustments](#rocket-heritability-adjustments)
* [Genetic correlation adjustments](#rocket-genetic-correlation-adjustments)
* [An example](#key-an-example)


## :white_check_mark: Prerequisites

The software is developed and tested in Linux and Windows environments.
- Python 2.7
- R (>=3.6)
- GNU Scientific Library (GSL) (>=2.3)

## :hammer_and_wrench: Installation
Download the jackknife LDSC software:

In R:
```r
devtools::install_github("shuangsong0110/ParticipationBias")
```

## :scroll: Prepare GWAS summary statistics
Please prepare the GWAS summary statistics in the following format (including the header line):
```
   chr        rsid     ref   alt       z         
    1      rs4040617    G     A     -0.199    
    1      rs4075116    C     T      0.646     
    1      rs9442385    T     G     -0.016    
    ...
```
**chr**: chromosome

**rsid**: SNP rsid

**ref**: reference allele

**alt**: alternative allele

**z**: GWAS z score



## :rocket: Heritability adjustments
### Step 1: Run LDSC


### Step 2: Making adjustments



## :rocket: Genetic correlation adjustments
### Step 1: Run LDSC


### Step 2: Making adjustments



## :key: An example
Download GWAS summary statistics for educational attainment:
```
wget
```

Perform LDSC:
```
ldsc.py
```

Make adjustments:
```
library(ParticipationBias)
```



