library(ggplot2)
library(dplyr)
library(purrr)
library(gridExtra)

res <- unlist(snakemake@input[["file"]])
prefix <- snakemake@params[["prefix"]]
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
out1 <- snakemake@output[["out1"]]
out2 <- snakemake@output[["out2"]]

res <- Filter(function(f) !is.null(readRDS(f)), res)
selection_methods <- res %>% strsplit(prefix) %>% sapply(tail, 1) %>%
  strsplit("_MVMR_") %>% sapply(head, 1)
selection_methods <- gsub('selection_', '', selection_methods) %>%
  strsplit("_seed") %>% sapply(head,1)
selection_methods[selection_methods=="unique_traits_filter"] <- "all"

all_res <- data.frame()
for (i in 1:length(res)) {
  sub_res <- readRDS(res[i])
  sub_res$selection_method <- selection_methods[i]
  all_res <- bind_rows(all_res,sub_res)
}


all_res <- all_res %>% mutate(CI_lower=b-qnorm(0.975)*se, CI_higher=b + qnorm(0.975)*se) %>%
           mutate(odds=exp(b),CI_lower=exp(CI_lower),CI_higher=exp(CI_higher))

# plot by odds
plt<- all_res %>% filter(exposure==id_exposure) %>%
       filter(converge == TRUE | is.na(converge)) %>%
       filter(se < 2) %>%
       ggplot() +
       geom_vline(xintercept = 1) +
       geom_point(aes(y = selection_method, x = odds, color = method,  group = method),
                  position=position_dodge(width = 0.9), size = 3) +
       geom_errorbar(aes(y = selection_method, xmin =CI_lower, xmax = CI_higher, color = method),
                     position=position_dodge(width = 0.9)) +
       xlab("Odds Ratio (95% CI)") + coord_flip() +
       theme_bw() + ggtitle(paste0("Direct causal effect of ",id_exposure,"  on ",id_outcome))+
       theme(axis.text.y = element_text(size = 20),
             axis.text.x = element_text(size = 15, angle = 10),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        legend.position = "bottom")
ggsave(out2, plot = plt,width = 20, height = 10)
write.csv(all_res,file = out1,row.names=FALSE)
