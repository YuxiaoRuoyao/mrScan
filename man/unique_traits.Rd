% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unique_traits.R
\name{unique_traits}
\alias{unique_traits}
\title{Calculate pairwise correlation and select unique traits}
\usage{
unique_traits(
  id.list,
  df_info,
  R_matrix,
  df_pair,
  R2_cutoff = 0.9,
  method = "cluster"
)
}
\arguments{
\item{id.list}{GWAS ID list of traits based on previous steps}

\item{df_info}{Dataframe of trait info from previous steps}

\item{R_matrix}{Pairwise correlation matrix}

\item{df_pair}{A dataframe contain string correlation for each pair}

\item{R2_cutoff}{high correlation cutoff to assign as duplicated traits. Default=.9}

\item{method}{filtering duplicate method: sample_size, nsnp, cluster. Default = "cluster"}
}
\value{
A GWAS ID vector and a trait info dataframe
}
\description{
Calculate pairwise correlation and select unique traits
}
