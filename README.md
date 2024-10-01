# ParticipationBias

A pipeline for adjusting participation bias in the estimation of heritability and genetic correlation.


## Table of contents
* [Prerequisites](#white_check_mark-prerequisites)
* [Installation](#hammer_and_wrench-installation)
* [Prepare GWAS summary statistics](#scroll-prepare-gwas-summary-statistics)
* [Example 1: Heritability adjustments](#rocket-example-1-heritability-adjustments)
* [Example 2: Genetic correlation adjustments](#rocket-example-2-genetic-correlation-adjustments)



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

## :rocket: Example 1: Heritability adjustments

Here we use the **heritability** of Educational Attainment (**EA**) as an example.
### Step 0: Download files
Download munged GWAS summary statistics for participation:
```
path=/home/local/ (change to your working path)
mkdir ./sumstats
wget -O ./sumstats/PB.sumstats.gz https://github.com/shuangsong0110/ParticipationBias/raw/refs/heads/main/example_data/PB.sumstats.gz
```

Download munged GWAS summary statistics for EA (Users could also specify their own GWAS summary statistics):
```
wget -O ./sumstats/EA.sumstats.gz https://github.com/shuangsong0110/ParticipationBias/raw/refs/heads/main/example_data/EA.sumstats.gz
```
### Step 1: Run LDSC (python 2)

```
path=/home/local/ (specify your working path)
trait_name='EA'
mkdir ${path}/results_${trait_name}
cd ${path}/results_${trait_name}
python ${path}/ldsc_jackknife/ldsc.py --rg ${path}/sumstats/PB.sumstats.gz,${path}/sumstats/${trait_name}.sumstats.gz --ref-ld ${path}/ldsc_jackknife/UKBB.EUR --w-ld ${path}/ldsc_jackknife/UKBB.EUR --intercept-gencov 0,0 --out res_rg
```

### Step 2: Making adjustments
In R:
```
library(ParticipationBias)
res_h2 <- h2_PB_adjust(path = '/home/local/', ## specify your working path, consistent to the LDSC path
                       mean_shift = 0.438,
                       trait_name = 'EA')
print(res_h2)
```


**path**: working path

**mean_shift**: mean shift of the phenotype of interest, between the sample of participant (UKBB) and the population, standardized in the sample of participants ((mean_participants-mean_population)/SE_in_participants)

**trait_name**: the name of the phenotype of interest

**trait_binary**: whether the trait is binary

**K**: prevalence of the binary trait



## :rocket: Example 2: Genetic correlation adjustments

Here we use the **genetic correlation** between Educational Attainment (**EA**) and **BMI** as an example.

### Step 0: Download files
Download munged GWAS summary statistics for **participation**:
```
path=/home/local/ (change to your working path)
mkdir ./sumstats
wget -O ./sumstats/PB.sumstats.gz https://github.com/shuangsong0110/ParticipationBias/raw/refs/heads/main/example_data/PB.sumstats.gz
```

Download munged GWAS summary statistics for **EA** (Users could also specify their own GWAS summary statistics):
```
wget -O ./sumstats/EA.sumstats.gz https://github.com/shuangsong0110/ParticipationBias/raw/refs/heads/main/example_data/EA.sumstats.gz
```

Download munged GWAS summary statistics for **BMI** (Users could also specify their own GWAS summary statistics):
```
wget -O ./sumstats/BMI.sumstats.gz https://github.com/shuangsong0110/ParticipationBias/raw/refs/heads/main/example_data/BMI.sumstats.gz
```

### Step 1: Run LDSC
**a. Participation & EA**
```
path=/home/local/ (specify your working path)
trait_name='EA'
mkdir ${path}/results_${trait_name}
cd ${path}/results_${trait_name}
python ${path}/ldsc_jackknife/ldsc.py --rg ${path}/sumstats/PB.sumstats.gz,${path}/sumstats/${trait_name}.sumstats.gz --ref-ld ${path}/ldsc_jackknife/UKBB.EUR --w-ld ${path}/ldsc_jackknife/UKBB.EUR --intercept-gencov 0,0 --out res_rg
```

**b. Participation & BMI**
```
path=/home/local/ (specify your working path)
trait_name='BMI'
mkdir ${path}/results_${trait_name}
cd ${path}/results_${trait_name}
python ${path}/ldsc_jackknife/ldsc.py --rg ${path}/sumstats/PB.sumstats.gz,${path}/sumstats/${trait_name}.sumstats.gz --ref-ld ${path}/ldsc_jackknife/UKBB.EUR --w-ld ${path}/ldsc_jackknife/UKBB.EUR --intercept-gencov 0,0 --out res_rg
```

**c. EA & BMI**
```
path=/home/local/ (specify your working path)
trait_name1='EA'
trait_name2='BMI'
mkdir ${path}/results_${trait_name1}_${trait_name2}
cd ${path}/results_${trait_name1}_${trait_name2}
python ${path}/ldsc_jackknife/ldsc.py --rg ${path}/sumstats/${trait_name1}.sumstats.gz,${path}/sumstats/${trait_name2}.sumstats.gz --ref-ld ${path}/ldsc_jackknife/UKBB.EUR --w-ld ${path}/ldsc_jackknife/UKBB.EUR --out res_rg
```

### Step 2: Making adjustments
In R:
```
library(ParticipationBias)
res_gcor <- gcor_PB_adjust(path = '/home/local/', ## specify your working path, consistent to the LDSC path
                      mean_shift1 = 0.438, mean_shift2 = -0.138,
                      trait_name1 = 'EA', trait_name2 = 'BMI')
print(res_gcor)
```

## :busts_in_silhouette: Maintainer

Please contact Shuang Song (shuangsong@hsph.harvard.edu) if there are any problems or questions.

## Acknowledgements
The GWAS summary statistics for BMI and EA are based on UKBB European samples.

The GWAS summary statistics for participation are derived with the method described in Benonisdottir and Kong (2023) (GWAS catalog accession codes: GCST90267220, GCST90267221, GCST90267222 and GCST90267223)

The original estimation of heritability and genetic correlation is based on LDSC method (https://github.com/bulik/ldsc/wiki).



