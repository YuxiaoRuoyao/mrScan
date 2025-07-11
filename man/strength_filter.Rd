% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strength_filter.R
\name{strength_filter}
\alias{strength_filter}
\title{Calculates the conditional F-statistic to assess instrument strength and filter traits}
\usage{
strength_filter(
  dat,
  dat_type = "local",
  R_matrix = NULL,
  df_info,
  pval_threshold = 5e-08,
  F_threshold = 5,
  effect_size_cutoff = 0.1,
  min_instruments = 3,
  Filter = FALSE,
  extra_traits = "None",
  type_outcome = "continuous",
  prevalence_outcome = NULL,
  type_exposure = NULL,
  prevalence_exposure = NULL,
  ss.exposure = NULL,
  df_af_out = NULL,
  df_af_exp = NULL
)
}
\arguments{
\item{dat}{A data frame of combined GWAS summary data after LD pruning.
For local data, the required columns include `SNP` for rsID, `trait_ID.beta` for beta hats,
`trait_ID.se` for standard errors, `trait_ID.p` for pvalues.
For IEU data, it should be output from TwoSampleMR::mv_harmonise_data().}

\item{dat_type}{Either "local" or "IEU". Default  = "local"}

\item{R_matrix}{Pairwise sample overlap matrix among traits. Default = NULL}

\item{df_info}{Dataframe of trait info from previous step}

\item{pval_threshold}{pvalue cutoff for selecting instruments. Default = 5e-8}

\item{F_threshold}{F-statistic cutoff. Default = 5}

\item{effect_size_cutoff}{Standardized effect size threshold. Default = 0.1}

\item{min_instruments}{minimum number of instruments. Default = 3}

\item{Filter}{Whether perform trait filtering based on F-stats and F_threshold. Default = FALSE}

\item{extra_traits}{trait ID you want to include no matter the instrument strength. Default = "None"}

\item{type_outcome}{It could be either "continuous" or "binary". Default = "continuous"}

\item{prevalence_outcome}{Outcome prevalence. It should be input if the outcome is a binary trait. Default = NULL}

\item{type_exposure}{A vector for the type of exposures. The order should be exactly matched
with exposures. eg. c("continuous","binary","continuous") for the second exposure is a binary trait}

\item{prevalence_exposure}{A vector for prevalence of exposures. The order should
be exactly matched with exposures. For continuous trait, just write NA. eg. c(NA, 0.1, NA)}

\item{ss.exposure}{A vector of sample size for exposures. You can provide it when dat_type = "IEU".
The order of it should be the same with beta hat matrix and se matrix. Default = NULL}

\item{df_af_out}{A dataframe of allele frequency of the outcome.
It contains columns `SNP`,`eaf.outcome`,`beta.outcome`,`id.outcome`. Default = NULL}

\item{df_af_exp}{A list for allele frequency matrix for each exposure.
Each dataframe contains columns `SNP`,`eaf.exposure`,`beta.exposure`,`id.exposure`.
Each element in the list is the dataframe for each exposure. Default = NULL}
}
\value{
A list of selected traits, a dataframe of conditional instrument strength and a dataframe of trait info
}
\description{
Calculates the conditional F-statistic to assess instrument strength and filter traits
}
