library(data.table)

current_id_list <- read.csv(snakemake@input[["new_id_list"]])$id
iteration_info_path <- snakemake@input[["iteration_file"]]

if (file.exists(iteration_info_path)) {
  current_iter <- read.csv(iteration_info_path)$iteration
} else {
  current_iter <- 0
  write.csv(data.frame(iteration = current_iter), iteration_info_path,row.name=FALSE)
}

prev_id_list <- if (current_iter > 0) {
    prev_iter <- current_iter - 1
    read.csv(snakemake@input[["prev_id_list"]])$id
} else {
    NULL  # For the very first iteration
}

if (!is.null(prev_id_list) && !identical(sort(current_id_list), sort(prev_id_list))) {
    write.csv(data.frame(continue = TRUE), snakemake@output[["continue_flag"]], row.names = FALSE, quote = FALSE)
    write.csv(data.frame(iteration = current_iter + 1), iteration_info_path, row.names = FALSE, quote = FALSE)
} else {
    file.copy(snakemake@input[["new_id_list"]], snakemake@output[["final_id_list"]], overwrite = TRUE)
}