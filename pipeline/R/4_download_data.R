library(stringr)
library(dplyr)
res <- readRDS(snakemake@input[["file"]])
df_harmonise <- read.csv(snakemake@input[["df_harmonise"]],header = F)
data_path <- snakemake@params[["path"]]
path_checkpoint <- snakemake@params[["checkpoint"]]
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
out <- snakemake@output[["out"]]

id_list <- unique(c(res$id.list,id_exposure,id_outcome))
ebi_list <- id_list[grep("GCST",id_list)]
regular_list <- id_list[!id_list %in% ebi_list]

GCST_list <- ebi_list %>% strsplit("-") %>% sapply(tail,1)
df_harmonise$V2 <- df_harmonise$V1 %>% strsplit("-") %>% sapply( "[", 3)
GCST_file <- df_harmonise %>% filter(str_detect(V2,paste0(GCST_list,collapse = '|'))) %>%
  pull(V1) %>% str_split("/",n=2) %>% sapply(tail,1)
nf_list <- ebi_list[!GCST_list %in% df_harmonise$V2]
if(length(nf_list) > 0){
  cat(paste0("Cannot find correct download link of ",nf_list,". Please input manually!"))
}
f1 <- paste0("wget -N -P ",data_path," https://gwas.mrcieu.ac.uk/files/",regular_list,"/",regular_list,".vcf.gz")
f2 <- paste0("wget -N -P ",data_path," https://gwas.mrcieu.ac.uk/files/",regular_list,"/",regular_list,".vcf.gz.tbi")
f3 <- paste0("wget -N -P ",data_path," https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",GCST_file)
checkpoint <- paste0("echo 'all done!' > ",path_checkpoint)
f <- data.frame(c(f1,f2,f3,checkpoint))
write.table(f,file = out,row.names = FALSE,
            col.names = FALSE, quote = FALSE)

