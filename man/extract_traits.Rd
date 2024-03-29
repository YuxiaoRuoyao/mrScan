% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_traits.R
\name{extract_traits}
\alias{extract_traits}
\title{Initially extract candidate traits}
\usage{
extract_traits(
  id_exposure,
  pval_x = 5e-08,
  pval_z = 1e-05,
  pop = "EUR",
  batch = c("ieu-a", "ieu-b", "ukb-b"),
  r2 = 0.001,
  kb = 10000,
  access_token = ieugwasr::check_access_token(),
  min_snps = 5,
  type_exposure = "IEU",
  type_candidate_traits = "IEU",
  file_path = NA,
  ref_path = NA,
  file_list = NA,
  trait_list = NA,
  snp_name_list = NA,
  beta_hat_name_list = NA,
  se_name_list = NA,
  p_value_name_list = NA
)
}
\arguments{
\item{id_exposure}{GWAS ID for main exposure.}

\item{pval_x}{p-value threshold to extract instruments of the main exposure. Default=5e-8}

\item{pval_z}{p-value threshold to retrieve traits. Default=1e-5}

\item{pop}{Super-population to use as reference panel. Default = "EUR". Options are "EUR","SAS","EAS","AFR","AMR"}

\item{batch}{GWAS database sub-batches, a vector. Default = c("ieu-a", "ieu-b","ukb-b")}

\item{r2}{clumping r2 threshold. Default=0.001}

\item{kb}{clumping kb window. Default=10000}

\item{access_token}{Google OAuth2 access token. Default=check_access_token()}

\item{min_snps}{the number of minimum shared SNPs with IV of X. Default=5}

\item{type_exposure}{Exposure data type. Either could be "IEU" or "local". Default = "IEU"}

\item{type_candidate_traits}{Candidate traits data type. Either could be "IEU" or "local". Default = "IEU"}

\item{file_path}{File path of local exposure GWAS summary data. It should be entered when you use local exposure data. Default = NA}

\item{ref_path}{LD reference data path. It should be entered when you use local exposure data. Default = NA}

\item{file_list}{GWAS summary data file paths for candidate traits. It should be entered when you use local traits data. Default = NA}
}
\value{
A GWAS ID vector and a trait info dataframe
}
\description{
Initially extract candidate traits
}
