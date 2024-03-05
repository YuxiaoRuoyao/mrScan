library(dplyr)
# cluster 49 ebi-a-GCST90029027, ebi-a-GCST90029026, ukb-b-7953
# cluster 95 ebi-a-GCST90013993, ebi-a-GCST90013994, ebi-a-GCST90025955
# cluster 254 ieu-b-4812, delete ieu-b-4808
# all cluster 111, 139, 177, 178, 179, 180, 181, 182, 189, 190, 193, 194, 195, 196, 207, 208, 221, 248, 262, 275, 278, 299
df_string <- read.csv("results/CRP_qc_string_trait_info.csv")
df_string_list <- read.csv("results/CRP_qc_string_id_list.csv")$id
id_list <- c("ebi-a-GCST90029027", "ebi-a-GCST90029026", "ukb-b-7953",
             "ebi-a-GCST90013993", "ebi-a-GCST90013994", "ebi-a-GCST90025955",
             "ieu-b-4812")
select_id <- df_string %>% filter(string_cluster %in%
                                    c(111, 139, 177, 178, 179, 180, 181, 182,
                                      189, 190, 193, 194, 195, 196, 207, 208,
                                      221, 248, 262, 275, 278, 299)) %>%
  pull(id) %>% vctrs::vec_c(id_list) %>% unique()
df_string <- df_string %>% filter(id %in% select_id) %>%
  mutate(status=="select after string similarity filtering")
df_string[df_string$id == "ieu-b-4808","status"] <- "delete due to string similarity"
df_string_list <- unique(c(df_string_list,select_id))
df_string_list <- df_string_list[df_string_list != "ieu-b-4808"]
write.csv(data.frame(id = df_string_list),file = "results/CRP_qc_string_id_list.csv",row.names = F)
write.csv(df_string,file = "results/CRP_qc_string_trait_info.csv",row.names = F)
