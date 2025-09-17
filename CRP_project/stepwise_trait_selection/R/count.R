library(dplyr)
library(stringr)

numbers <- read.csv(snakemake@input[["id_list"]])$id
out_id_list <- snakemake@output[["out_id_list"]]
out_flag <- snakemake@output[["out_flag"]]


if (length(numbers) > 1) {
  numbers <- numbers[numbers != min(numbers)]
  write.csv(data.frame(id = numbers), out_id_list,row.names=F)
  write("CONTINUE", out_flag)
} else {
  write.csv(data.frame(id = numbers), out_id_list,row.names=F)
  write("STOP", out_flag)
}