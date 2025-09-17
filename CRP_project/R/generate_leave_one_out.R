
id_list <- read.csv(snakemake@input[["id_list"]])$id
leaveout_id <- snakemake@params[["leaveout_id"]]
out <- snakemake@output[["out"]]

keep_ids <- id_list[id_list != leaveout_id]
write.csv(data.frame(id = keep_ids), file = out, row.names=FALSE)