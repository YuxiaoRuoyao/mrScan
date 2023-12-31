% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_R_pval.R
\name{estimate_R_pval}
\alias{estimate_R_pval}
\title{Estimate genetic correlation matrix by pvalue}
\usage{
estimate_R_pval(
  ld_prune_file_dir = NULL,
  prefix = NULL,
  p_thresh = 0.05,
  out_dir = NULL
)
}
\arguments{
\item{prefix}{Name prefix for the output. Default = NULL}

\item{out_dir}{Output data path. Default is in the current work directory.}

\item{beta_dir}{Directory path of merged data. Default is the current work directory.}

\item{ref_path}{Path for the LD reference panel.}
}
\value{
Save one dataframe per chromosome with columns for SNP info
}
\description{
Estimate genetic correlation matrix by pvalue
}
