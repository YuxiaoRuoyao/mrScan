library(dplyr)
# cluster 95 ebi-a-GCST90013993, ebi-a-GCST90013994, ebi-a-GCST90025955
# all cluster 139, 246, 260, 273, 296, 299
# cluster 235 only ebi-a-GCST90104006
# cluster 252 ieu-b-4812
df_string <- read.csv("results/CRP_qc_string_trait_info.csv")
df_string_list <- read.csv("results/CRP_qc_string_id_list.csv")$id
id_list <- c("ebi-a-GCST90013993", "ebi-a-GCST90013994", "ebi-a-GCST90025955","ieu-b-4812")
select_id <- df_string %>% filter(string_cluster %in% c(139, 246, 260, 273, 296, 299)) %>%
  pull(id) %>% vctrs::vec_c(id_list) %>% unique()
df_string[df_string$id %in% select_id, "status"] <- "select after string similarity filtering"
df_string[df_string$id == "ieu-b-4808","status"] <- "delete due to string similarity"
df_string_list <- unique(c(df_string_list,select_id))
df_string_list <- df_string_list[df_string_list != "ieu-b-4808"]
id_235 <- df_string %>% filter(string_cluster == 235 & id != "ebi-a-GCST90104006") %>% pull(id)
df_string_list <- df_string_list[!df_string_list %in% id_235]
df_string[df_string$id %in% id_235,"status"] <- "delete due to string similarity"

write.csv(data.frame(id = df_string_list),file = "results/CRP_qc_string_id_list.csv",row.names = F)
write.csv(df_string,file = "results/CRP_qc_string_trait_info.csv",row.names = F)
