% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confounder_selection.R
\name{confounder_selection}
\alias{confounder_selection}
\title{Select confounders by variable selection}
\usage{
confounder_selection(
  id_exposure,
  id_outcome,
  id.list,
  df_info,
  method,
  stepwise_method = "forward",
  r2 = 0.001,
  kb = 10000,
  pval_threshold = 5e-08,
  find_proxies = TRUE,
  pop = "EUR",
  harmonise_strictness = 2,
  lambda_type = "min",
  radius_type = "1se",
  seed = 1,
  maxits = 1e+06,
  p_cutoff = 0.05,
  df_bidirection = NULL
)
}
\arguments{
\item{id_exposure}{GWAS ID of the main exposure}

\item{id_outcome}{GWAS ID of the outcome}

\item{id.list}{GWAS ID list of traits based on downstream_filter}

\item{df_info}{Dataframe of trait info from previous steps}

\item{method}{Variable selection method. Include Lasso-type, marginal selection and self-input}

\item{stepwise_method}{either forward, backward or both. Default = "forward"}

\item{r2}{LD-clump r2. Default = 0.001}

\item{kb}{LD-clump kb. Default = 10000}

\item{pval_threshold}{Instrument selection cutoff. Default = 5e-8}

\item{find_proxies}{Whether look for proxies. Default = TRUE}

\item{pop}{Super population to use. Default = "EUR"}

\item{harmonise_strictness}{Data harmonise strictness. See documentation of TwoSample MR}

\item{lambda_type}{Lasso penalty type. Default = "min"}

\item{radius_type}{Radius type for corresponding loss to corrected Lasso. Default = "min"}

\item{seed}{Default = 1}

\item{maxits}{maximum number of iterations for algorithm converges. Default=1000000}

\item{p_cutoff}{pvalue threshold for marginal selection}

\item{df_bidirection}{dataframe of traits with estimates of four direction, used
by marginal selection. Results from downstream_filter}
}
\value{
A GWAS ID vector and a trait info dataframe
}
\description{
Select confounders by variable selection
}
