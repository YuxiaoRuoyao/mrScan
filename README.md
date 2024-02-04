
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mrScan

This is the R package for mrScan (Automatically Select Heritable
Confounders for Mendelian Randomization). This package will help you to
find phenome-wide potential heritable confounders and give direct causal
estimates of the main exposure to the outcome after adjusting for
confounders by multivariable Mendelian Randomization (MVMR).

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
devtools::install_github("jean997/esmr")
```

Please note the vignette only provides a simple example of function
usage. If you want to do a systematically searching for heritable
confounders and accurately get the causal estimate of the main exposure,
we recommend you to use the prepared Snakemake pipeline.

## Snakemake pipeline usage

### Install Snakemake

See the official install instruction
[here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
You can also quickly learn the usage of Snakemake
[here](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

### Install necessary R packages

You probably will need the following R packages if you want to use local
GWAS summary statistics.

``` r
devtools::install_github("mrcieu/gwasvcf")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
install.packages("rlang")
remotes::install_github("privefl/bigsnpr")
```

### Download the Snakemake pipeline

You need to download the whole pipeline directory to your local space.
There are a few ways to that:

- Use [Download Directory](https://download-directory.github.io/). Just
  copy address link
  <https://github.com/YuxiaoRuoyao/mrScan/tree/master/pipeline> and
  download.  
- Manually use command line via SVN.

``` r
brew install svn
svn checkout https://github.com/YuxiaoRuoyao/mrScan/trunk/pipeline
```

### (Optional Step) Create a csv file with all candidate traits if you want to use local data

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
