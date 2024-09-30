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
wget -O ldsc_jackknife.tar.gz "https://dl.dropboxusercontent.com/scl/fi/3lgslbgqz4c1sebje0473/ldsc_jackknife.tar.gz?rlkey=5l6c0mwgljamnbs3ddpearu2v&st=26gjcp8s&dl=1"
tar -zxvf ldsc_jackknife.tar.gz
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
path=/home/local/ (change to your working path)
mkdir ./sumstats
wget -O ./sumstats/PB.sumstats.gz https://github.com/shuangsong0110/ParticipationBias/raw/refs/heads/main/example_data/PB.sumstats.gz
```

Download munged GWAS summary statistics for educational attainment:
```
wget -O ./sumstats/EA.sumstats.gz https://github.com/shuangsong0110/ParticipationBias/raw/refs/heads/main/example_data/EA.sumstats.gz
```

Perform LDSC (python2):
```
path=/home/local/ (specify your working path)
trait_name='EA'
mkdir ./results_${trait_name}
cd ./results_${trait_name}
python ${path}/ldsc_jackknife/ldsc.py --rg ${path}/sumstats/PB.sumstats.gz,${path}/sumstats/${trait_name}.sumstats.gz --ref-ld ${path}/ldsc_jackknife/UKBB.EUR --w-ld ${path}/ldsc_jackknife/UKBB.EUR --intercept-gencov 0,0 --out res_rg
```

Make adjustments:
```
library(ParticipationBias)
res <- h2_PB_adjust(path='/home/local/', ## specify your working path, consistent to the LDSC path
                    mean_shift=0.438,
                    trait_name='EA')
print(res)
```



