% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downstream_filter.R
\name{downstream_filter}
\alias{downstream_filter}
\title{Filter downstream traits}
\usage{
downstream_filter(
  id_exposure,
  id_outcome,
  id.list,
  df_info,
  method = c("mr_raps"),
  sig_level = 0.05
)
}
\arguments{
\item{id_exposure}{GWAS ID of the main exposure}

\item{id_outcome}{GWAS ID of the outcome}

\item{id.list}{GWAS ID list of traits based on previous steps}

\item{df_info}{Dataframe of trait info from previous steps}

\item{method}{Bidirection MR method. Same with TwoSampleMR package.Default = "mr_raps"}

\item{sig_level}{One-sided t-test significant level. Default=0.05}
}
\value{
A GWAS ID vector, a trait info dataframe, a trait dataframe with four direction
estimate and t-test results
}
\description{
Filter downstream traits
}
