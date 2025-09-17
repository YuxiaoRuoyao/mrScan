id_list <- read.csv(snakemake@input[["id_list"]])$id
trait_info <- read.csv(snakemake@input[["trait_info"]])
out <- snakemake@output[["out"]]

trait_info[trait_info$id %in% id_list,"status"] <- "select after stepwise strength filtering"
write.csv(trait_info, file = out, row.names = F)