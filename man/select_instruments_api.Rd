% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_instruments_api.R
\name{select_instruments_api}
\alias{select_instruments_api}
\title{Select instruments by api for MVMR analysis}
\usage{
select_instruments_api(
  id.list,
  id_exposure,
  id_outcome,
  r2 = 0.001,
  kb = 10000,
  pval_threshold = 5e-08,
  find_proxies = TRUE,
  pop = "EUR",
  harmonise_strictness = 2
)
}
\arguments{
\item{id.list}{GWAS ID list of traits based on previous steps}

\item{id_exposure}{GWAS ID of the main exposure}

\item{id_outcome}{GWAS ID of the outcome}

\item{r2}{LD-clump r2. Default = 0.001}

\item{kb}{LD-clump kb. Default = 10000}

\item{pval_threshold}{Instrument selection cutoff. Default = 5e-8}

\item{find_proxies}{Whether look for proxies. Default = TRUE}

\item{pop}{Super population to use. Default = "EUR"}

\item{harmonise_strictness}{Data harmonise strictness. See documentation of TwoSample MR}

\item{df_info}{Dataframe of trait info from previous steps}
}
\value{
A list of harmonized data for the outcome (mvdat_y) and a list of harmonized data for the main exposure (mvdat_x)
}
\description{
Select instruments by api for MVMR analysis
}
