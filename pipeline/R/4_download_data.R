res <- readRDS(snakemake@input[["file"]])
data_path <- snakemake@params[["path"]]
out <- snakemake@output[["out"]]

id_list <- res$id.list
f1 <- paste0("wget -N -P ",data_path," https://gwas.mrcieu.ac.uk/files/",id_list,"/",id_list,".vcf.gz")
f2 <- paste0("wget -N -P ",data_path," https://gwas.mrcieu.ac.uk/files/",id_list,"/",id_list,".vcf.gz.tbi")
f <- data.frame(c(f1,f2))
write.table(f,file = out,row.names = FALSE,
            col.names = FALSE, quote = FALSE)

