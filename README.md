
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Installation with vignettes

You can install the development version of mrScan from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("YuxiaoRuoyao/mrScan")
#browseVignettes("mrScan")
```

You should also install the following packages in advance for vignette:

``` r
install.packages("dplyr")
install.packages("remotes")
install.packages("stringr")
remotes::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("mrcieu/ieugwasr")
devtools::install_github("qingyuanzhao/mr.raps")
devtools::install_github("noahlorinczcomi/MRBEE")
devtools::install_github("jean997/sumstatFactors")
devtools::install_github("jingshuw/grapple")
```

## (Optional Step) Create a csv file with all candidate traits if you want to use local data

If you have a list of potential candidate traits and have already
downloaded GWAS summary data to local, you need to create a csv file
containing the file information of local data. The GWAS data itself can
be in one of three formats:

- A vcf file downloaded from the IEU Open GWAS database
- A harmonized .h.tsv.gz file downloaded from EBI database
- A flat file with columns for snp, effect and non-effect alleles,
  effect estimate, and standard error

The csv file should contain the following columns and use NA to indicate
missing data:

- trait_ID: The GWAS trait ID for this specific study. Basically just
  use the trait ID in the IEU database.
- path: Local path to the raw data files.  
- snp: Column name for rsid.  
- beta_hat: Column name for coefficient estimate.
- se: Column name for standard error of beta_hat.
- p_value: Column name for p-value.

Note the pvalue column content could be NA if the data does not have
pvalue. If the local file is downloaded from IEU Open GWAS database or
EBI database (harmonized one), you can just write NA for snp, beta_hat,
se and p_value column and don’t need to check the exact column name. But
you must need to provide snp, beta_hat and se column names if it’s a
flat file because we don’t know your own data structure in advance.
