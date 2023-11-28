
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
#> Installing package into '/private/var/folders/09/lsxqw2292sld7jrjxdslj7_w0000gn/T/RtmpIvTxHM/temp_libpath152854be5cc70'
#> (as 'lib' is unspecified)
#> 
#>   There is a binary version available but the source version is later:
#>       binary source needs_compilation
#> dplyr  1.1.3  1.1.4              TRUE
#> installing the source package 'dplyr'
install.packages("remotes")
#> Installing package into '/private/var/folders/09/lsxqw2292sld7jrjxdslj7_w0000gn/T/RtmpIvTxHM/temp_libpath152854be5cc70'
#> (as 'lib' is unspecified)
#> 
#> The downloaded binary packages are in
#>  /var/folders/09/lsxqw2292sld7jrjxdslj7_w0000gn/T//Rtmp455CIz/downloaded_packages
remotes::install_github("MRCIEU/TwoSampleMR")
#> Skipping install of 'TwoSampleMR' from a github remote, the SHA1 (cbd03e6a) has not changed since last install.
#>   Use `force = TRUE` to force installation
devtools::install_github("qingyuanzhao/mr.raps")
#> Skipping install of 'mr.raps' from a github remote, the SHA1 (2a23d841) has not changed since last install.
#>   Use `force = TRUE` to force installation
```
