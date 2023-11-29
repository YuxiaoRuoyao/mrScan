#' @title Format GWAS summary statistics from ieu-OpenGWAS database
#' @param data_file Raw GWAS summary data path
#' @param output_file Output data path
#' @returns Formatted GWAS summary data with .vcf.bgz
#'
#' @import stringr
#' @import readr
#' @import VariantAnnotation
#' @import gwasvcf
#' @import dplyr
#' @import magrittr
#' @export
format_gwas_data <- function(data_file,output_file){
    #Harmonize strand in vcf to match cause format (A1 = A)
    # Note REF = A2, ALT = A1
    dat <- readVcf(data_file) %>%
      vcf_to_tibble()
    dat <- dat %>%
      rename(A1 = ALT, A2 = REF)
    #remove non snp
    l2 <- str_length(dat$A2)
    l1 <- str_length(dat$A1)
    dat <- dat[l1==1 & l2==1,]
    dat <- sumstatFactors:::remove_ambiguous(dat)
    dat1 <- sumstatFactors:::align_beta(dat, "ES")

    dat <- dat1 %>%
      mutate(AF = case_when(ES == -1*dat$ES ~ 1-AF,
                            TRUE ~ AF)) %>%
      rename(ALT = A1, REF = A2) %>%
      filter(!is.na(ID))
    out <- dat %$% create_vcf(chrom=seqnames, pos=start,nea=REF,
                              ea=ALT, snp=ID, ea_af=AF,
                              effect=ES,  se=SE, pval=10^-LP, n=SS, name="a")
    output_file <- str_replace(output_file, ".bgz$", "")
    writeVcf(out, file=output_file, index=TRUE)
}


