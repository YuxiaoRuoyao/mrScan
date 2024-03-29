% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal.R
\name{marginal}
\alias{marginal}
\title{Use marginal selection to do confounder selection}
\usage{
marginal(
  id.list,
  df_info,
  df_bidirection,
  extra_traits = "None",
  p_cutoff = 0.05
)
}
\arguments{
\item{id.list}{GWAS ID list of traits based on previous steps}

\item{df_info}{Dataframe of trait info from previous steps}

\item{df_bidirection}{Dataframe of trait with four direction estimates. Result from downstream_filter}

\item{extra_traits}{trait ID which is adjusted for in bidirection MR. Default = "None"}

\item{p_cutoff}{pvalue threshold. Default = 0.05}
}
\value{
A GWAS ID vector and a trait info dataframe
}
\description{
Use marginal selection to do confounder selection
}
