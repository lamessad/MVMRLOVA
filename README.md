
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MVMRLOVA

<!-- badges: start -->
<!-- badges: end -->

The goal of MVMRLOVA is to performs multivariable Mendelian
randomization analysis based on latent phenotype of the outcome that
explicitly exclude vertical pleiotropy effects via EM algorithm. There
are three functions:

- **`MVMRLOVA()`** main function that performs MVMR analysis and
  provides causal effect estimates.

## Installation

You can install the development version of **MVMRLOVA** from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lamessad/MVMRLOVA")
```

To view the help pages for functions in this package, prepend the
function name with a question mark.

``` r
#library(MVMRLOVA)
#?MVMRLOVA()
```

## Example

Here is an example which shows how to apply MVMR-cML to infer the causal
relationship from cholesterol levels (triglycerides, LDL and HDL) to
CAD.

First extract GWAS summary data with package:

``` r
#IVs extraction
library(TwoSampleMR)
#> TwoSampleMR version 0.6.9
#>   [>] New authentication requirements: https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication.
#>   [>] Major upgrades to our servers completed to improve service and stability.
#>   [>] We need your help to shape our emerging roadmap!
#>       Please take 2 minutes to give us feedback -
#>       https://forms.office.com/e/eSr7EFAfCG
#> 
#> Warning:
#> You are running an old version of the TwoSampleMR package.
#> This version:   0.6.9
#> Latest version: 0.6.18
#> Please consider updating using remotes::install_github('MRCIEU/TwoSampleMR')
exposuregwas_id = c(
  'ebi-a-GCST002216',
  'ebi-a-GCST002222')
outcomegwas_id = 'ebi-a-GCST005195'
exposure_dat = mv_extract_exposures(exposuregwas_id)
#> Please look at vignettes for options on running this locally if you need to run many instances of this command.
#> Clumping 1, 126 variants, using EUR population reference
#> Removing 25 of 126 variants due to LD with other variants or absence from LD reference panel
#> Extracting data for 101 SNP(s) from 2 GWAS(s)
#> Harmonising Triglycerides || id:ebi-a-GCST002216 (ebi-a-GCST002216) and LDL cholesterol || id:ebi-a-GCST002222 (ebi-a-GCST002222)
#> Removing the following SNPs for being palindromic with intermediate allele frequencies:
#> rs3758348, rs7534572
outcome_dat <- extract_outcome_data(exposure_dat$SNP,outcomegwas_id)
#> Extracting data for 101 SNP(s) from 1 GWAS(s)
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
#> Harmonising Triglycerides || id:ebi-a-GCST002216 (ebi-a-GCST002216) and Coronary artery disease || id:ebi-a-GCST005195 (ebi-a-GCST005195)
#> Removing the following SNPs for being palindromic with intermediate allele frequencies:
#> rs3758348, rs7534572
betaX   <- mvdat$exposure_beta
colnames(betaX) <- c("TG", "LDL")
betaXse <- mvdat$exposure_se
betaY   <- as.matrix(mvdat$outcome_beta)  # make it a column matrix
betaYse <- as.matrix(mvdat$outcome_se)
ny=min(outcome_dat$samplesize.outcome)
#MVMRLOVA analysis
library(MVMRLOVA)
est <- MVMRLOVA(betaY, betaX, betaYse, betaXse, ny, gwas_p = 0.05, gwas_p2 = 5e-8, 100)
#> Warning: package 'MendelianRandomization' was built under R version 4.3.3
#> 
#> Attaching package: 'MendelianRandomization'
#> The following objects are masked from 'package:TwoSampleMR':
#> 
#>     mr_ivw, mr_median
#> Warning in MVMRLOVA(betaY, betaX, betaYse, betaXse, ny, gwas_p = 0.05, gwas_p2
#> = 5e-08, : To get a more precise p-value, it is recommended to increase the
#> number of permutations.
est$CausEst
#>      BxTG     BxLDL 
#> 0.2143080 0.2848379
```
