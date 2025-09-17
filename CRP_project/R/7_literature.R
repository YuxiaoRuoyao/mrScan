library(ieugwasr)
library(dplyr)

id.list <- read.csv(snakemake@input[["id_list"]])$id
df_info <- read.csv(snakemake@input[["trait_info"]])
literature_traits <- snakemake@params[["literature_traits"]]
out_id_list <- snakemake@output[["out_id_list"]]
out_trait_info <- snakemake@output[["out_trait_info"]]

id_literature <- unique(c(id.list,literature_traits))
new_info <- df_info
if(sum(!literature_traits %in% df_info$id) != 0){
  new_id <- literature_traits[!literature_traits %in% df_info$id]
  new_info <- df_info %>% full_join(gwasinfo(new_id))
}
new_info[new_info$id %in% id_literature,"status"] <- "Select by literature"

write.csv(data.frame(id = id_literature),file = out_id_list,row.names = F)
write.csv(new_info,file = out_trait_info,row.names = F)
