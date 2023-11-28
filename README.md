
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Installation with vignettes

You can install the development version of mrScan from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("YuxiaoRuoyao/mrScan",build_vignettes = TRUE)
browseVignettes("mrScan")
```

You should also install the following packages in advance for vignette:

``` r
install.packages("dplyr")
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("qingyuanzhao/mr.raps")
```
