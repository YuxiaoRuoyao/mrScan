% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MVMR_ESMR.R
\name{MVMR_ESMR}
\alias{MVMR_ESMR}
\title{Use ESMR to do MVMR analysis by locally data}
\usage{
MVMR_ESMR(
  dat,
  R_matrix,
  pval_threshold = 5e-08,
  effect_size_cutoff = 0.1,
  type_outcome = "continuous",
  prevalence_outcome = NULL,
  type_exposure = NULL,
  prevalence_exposure = NULL
)
}
\arguments{
\item{dat}{A data frame of combined GWAS summary data after LD pruning.
The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats,
`trait_ID.se` for standard errors, `trait_ID.z` for z-values, `trait_ID.ss` for sample sizes.}

\item{R_matrix}{Pairwise sample overlap matrix among traits}

\item{pval_threshold}{pvalue cutoff for selecting instruments. Default = 5e-8}

\item{effect_size_cutoff}{Standardized effect size threshold. Default = 0.1}

\item{type_outcome}{It could be either "continuous" or "binary". Default = "continuous"}

\item{prevalence_outcome}{Outcome prevalence. It should be input if the outcome is a binary trait. Default = NULL}

\item{type_exposure}{A vector for the type of exposures. The order should be exactly matched
with exposures. eg. c("continuous","binary","continuous") for the second exposure is a binary trait}

\item{prevalence_exposure}{A vector for prevalence of exposures. The order should
be exactly matched with exposures. For continuous trait, just write NA. eg. c(NA, 0.1, NA)}
}
\value{
A dataframe of result summary
}
\description{
Use ESMR to do MVMR analysis by locally data
}
