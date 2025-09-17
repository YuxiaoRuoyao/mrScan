plus_id <- snakemake@params[["plus_id"]]
out <- snakemake@output[["out"]]

keep_ids <- c("ebi-a-GCST90029070","ukb-b-19953",plus_id)
write.csv(data.frame(id = keep_ids), file = out, row.names=FALSE)