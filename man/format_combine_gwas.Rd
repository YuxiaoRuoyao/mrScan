% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_combine_gwas.R
\name{format_combine_gwas}
\alias{format_combine_gwas}
\title{Combine and harmonize local GWAS summary data with multiple traits}
\usage{
format_combine_gwas(df_file, c, df_info)
}
\arguments{
\item{df_file}{A dataframe contain "id" (trait id) and "location" (GWAS summary data location) columns.}

\item{c}{Numeric chromosome number}

\item{df_info}{Dataframe of trait info from previous steps containing sample size}
}
\value{
Save one dataframe per chromosome with columns for SNP info
}
\description{
Combine and harmonize local GWAS summary data with multiple traits
}
