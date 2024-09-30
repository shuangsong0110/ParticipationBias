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
```
wget -O
gzip -d 
```

In R:
```r
devtools::install_github("shuangsong0110/ParticipationBias")
```

## :scroll: Prepare GWAS summary statistics
Please prepare the GWAS summary statistics in the following format (including the header line, sep='\t'):
```
     SNP      A1    A2       Z         N       
 rs4040617    G     A     -0.199     360000
 rs4075116    C     T      0.646     360000
 rs9442385    T     G     -0.016     360000
    ...
```

**SNP**: SNP rsid

**A1**: reference allele

**A2**: alternative allele

**Z**: GWAS z score

**N**: GWAS sample size

### Munge summary statistics for trait 1:
```
python2 ./munge_sumstats.py --sumstats ./trait1.txt  --merge-alleles pan.snipar.snplist --out ./trait1.summs
```

## :rocket: Heritability adjustments
### Step 1: Run LDSC
```
cd ${result_path}
python ${ldsc_path}/ldsc.py --rg ${summstats_path}/PB.sumstats.gz,${summstats_path}/${trait_name}.sumstats.gz --ref-ld ${ldsc_path}/UKBB.EUR --w-ld ${ldsc_path}/UKBB.EUR --intercept-gencov 0,0
```

### Step 2: Making adjustments
```
library(ParticipationBias)
res <- h2_PB_adjust(path, mean_shift, trait_name='trait1', trait_binary=F, K=1)
```
**path**: working path

**mean_shift**: mean shift of the phenotype of interest, between the sample of participant (UKBB) and the population, standardized in the sample of participants ((mean_participants-mean_population)/SE_in_participants)

**trait_name**: the name of the phenotype of interest

**trait_binary**: whether the trait is binary

**K**: prevalence of the binary trait



## :rocket: Genetic correlation adjustments
### Step 1: Run LDSC
```
cd ${result_path}
python ${ldsc_path}/ldsc.py --rg ${summstats_path}/PB.sumstats.gz,${summstats_path}/${trait_name1}.sumstats.gz --ref-ld ${ldsc_path}/UKBB.EUR --w-ld ${ldsc_path}/UKBB.EUR --intercept-gencov 0,0
python ${ldsc_path}/ldsc.py --rg ${summstats_path}/PB.sumstats.gz,${summstats_path}/${trait_name2}.sumstats.gz --ref-ld ${ldsc_path}/UKBB.EUR --w-ld ${ldsc_path}/UKBB.EUR --intercept-gencov 0,0
${summstats_path}/${trait_name1}.sumstats.gz,${summstats_path}/${trait_name2}.sumstats.gz --ref-ld ${ldsc_path}/UKBB.EUR --w-ld ${ldsc_path}/UKBB.EUR 
```

### Step 2: Making adjustments



## :key: An example
Download munged GWAS summary statistics for participation:
```
wget -O PB.sumstats.gz https://hu-my.sharepoint.com/:u:/r/personal/shuangsong_hsph_harvard_edu/Documents/research_share/ParticipationBias/PB.sumstats.gz?csf=1&web=1&e=OtgQ1L
```

Download munged GWAS summary statistics for educational attainment:
```
wget -O EA.sumstats.gz https://hu-my.sharepoint.com/:u:/r/personal/shuangsong_hsph_harvard_edu/Documents/research_share/ParticipationBias/EA.sumstats.gz?csf=1&web=1&e=bXlEGz
```

Perform LDSC (python2):
```
python ./ldsc.py --rg PB.sumstats.gz,EA.sumstats.gz --ref-ld UKBB.EUR --w-ld UKBB.EUR --intercept-gencov 0,0
```

Make adjustments:
```
library(ParticipationBias)
```



