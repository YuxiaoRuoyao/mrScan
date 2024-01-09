library(ieugwasr)
library(dplyr)
library(data.table)
library(purrr)
mr_files <- unlist(snakemake@input[["mr_files"]])
id_file <- read.csv(snakemake@input[["id_list"]])
df_info <- read.csv(snakemake@input[["trait_info"]])
id_exposure <- snakemake@params[["id_exposure"]]
sig_level <- as.numeric(snakemake@params[["sig_level"]])
R2_cutoff <- as.numeric(snakemake@params[["R2_cutoff"]])
extra_trait <- snakemake@params[["extra_trait"]]
out <- snakemake@output[["out"]]

id.list <- id_file$id
res_mr <- map(mr_files, function(f){
  readRDS(f)
})
res <- do.call(Map, c(f = rbind, res_mr))

df_summary<-data.frame(id=id.list) %>% left_join(df_info[,c("id","trait")],by="id")
df_YtoZ <- left_join(df_summary,res$mr12[,c("id.exposure","id.outcome","b","se","pval")],
                     by=c("id"="id.outcome"))
df_YtoZ <- df_YtoZ %>% mutate(exposure = if_else(id.exposure == id_exposure, "X", "Y"))
setDT(df_YtoZ)
data_wide1 <- dcast(df_YtoZ,id ~ exposure, value.var=c("b","se")) %>%
              select_if(~sum(!is.na(.)) > 0)
colnames(data_wide1)<-c("id","b_XtoZ","b_YtoZ","se_XtoZ","se_YtoZ")
df_ZtoY <- left_join(df_summary,res$mr21[,c("id.exposure","id.outcome","b","se","pval")],
                     by=c("id"="id.exposure"))
df_ZtoY <- df_ZtoY %>% mutate(outcome = if_else(id.outcome == id_exposure, "X", "Y"))
setDT(df_ZtoY)
data_wide2 <- dcast(df_ZtoY,id ~ outcome, value.var=c("b","se")) %>%
              select_if(~sum(!is.na(.)) > 0)
colnames(data_wide2)<-c("id","b_ZtoX","b_ZtoY","se_ZtoX","se_ZtoY")
df_final<-left_join(data_wide1,data_wide2,by="id")
df_final<- df_final %>%
  mutate(t_Z_X = (abs(b_ZtoX)-abs(b_XtoZ))/sqrt(se_ZtoX^2+se_XtoZ^2),
         t_Z_Y = (abs(b_ZtoY)-abs(b_YtoZ))/sqrt(se_ZtoY^2+se_YtoZ^2)) %>%
  mutate(p_Z_X = pnorm(t_Z_X),p_Z_Y = pnorm(t_Z_Y))
X_downstream<- df_final %>% filter(p_Z_X < sig_level) %>% pull(id)
X_sig<- df_final %>% filter(p_Z_X > 1 - sig_level) %>% pull(id)
Y_downstream<- df_final %>% filter(p_Z_Y < sig_level) %>% pull(id)
Y_sig<- df_final %>% filter(p_Z_Y > 1 - sig_level) %>% pull(id)
sig_traits <- union(X_sig,Y_sig)
downstream_traits <- union(X_downstream,Y_downstream)
select_trait <- sig_traits[!sig_traits %in% downstream_traits]
df_info[df_info$id %in% select_trait,"status"] <- "select after downstream filtering"
filter.trait <- id.list[!id.list %in% select_trait]
df_info[df_info$id %in% filter.trait,"status"] <- "delete in downstream filtering"
# delete high correlation traits with either X and Y
res_cor <- res$cor
trait_cor_X <- res_cor %>% filter(id1 == id_exposure) %>%
  filter(abs(cor) > R2_cutoff) %>%
  pull(id2)
trait_cor_Y <- res_cor %>% filter(id1 != id_exposure) %>%
  filter(abs(cor) > R2_cutoff) %>%
  pull(id2)
df_info[df_info$id %in% trait_cor_X,"status"] <- "delete since high cor with X"
df_info[df_info$id %in% trait_cor_Y,"status"] <- "delete since high cor with Y"
select_trait <- select_trait[!select_trait %in% c(trait_cor_X,trait_cor_Y)]
if(extra_trait != "None"){
  select_trait <- c(select_trait,extra_trait)
  df_info[df_info$id %in% extra_trait,"status"] <- "select after downstream filtering"
}

saveRDS(list(id.list=select_trait,trait.info=df_info,df_bidirection = df_final),
        file = out)
