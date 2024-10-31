library(ieugwasr)
library(dplyr)

Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ5dXhpYW93QHVtaWNoLmVkdSIsImlhdCI6MTcyOTY0MDQ5NSwiZXhwIjoxNzMwODUwMDk1fQ.YZA2h58w7WgJjxWzNnKKnYRz84c7Gf40FVTP4aKMdYpa84DkcPk-IoHFeJgjMEeIhwGEC7gXnrtb6MTW0pp24jEESZFVX7_hIg5k0UnX3q6Y4o9GhjppNo5dZFof9ZRZmMzu43WCrUj_nYNcWL-k7vGaXhpGaS2AC1sAk24RLDHY23HYyPEQ76O0Cs4hBDjF2veNj8XAbV_6ThSLX8NA2g4ygd3FiZjlXk8cEPeEF49thNfNQWOJD3rcfI52V5AIn1VJk5gIN4PD-fm1GLWHou02EDdlVt_-sewZlPjmkpu7jVMjP9OZnppwIYj4w6Bi6uf4uYeFzMyiIrP6bKYPng")

id.list <- read.csv(snakemake@input[["id_list"]])$id
df_info <- read.csv(snakemake@input[["trait_info"]])
literature_traits <- snakemake@params[["literature_traits"]]
out_id_list <- snakemake@output[["out_id_list"]]
out_trait_info <- snakemake@output[["out_trait_info"]]

id_literature <- unique(c(id.list,literature_traits))
if(sum(!literature_traits %in% df_info$id) != 0){
  new_id <- literature_traits[!literature_traits %in% df_info$id]
  new_info <- df_info %>% full_join(gwasinfo(new_id))
}
new_info[new_info$id %in% id_literature,"status"] <- "Select by literature"

write.csv(data.frame(id = id_literature),file = out_id_list,row.names = F)
write.csv(new_info,file = out_trait_info,row.names = F)
