library(ggplot2)
library(dplyr)
all_res_0.05 <- data.frame()
all_res_0.1 <- data.frame()
all_res_MR <- data.frame()
# Generate all plot
label <- c("SCZ","T2D","CAD","stroke","Alzheimer","Knee","Colorectal","Parkinson",
           "IBD","BD","BMD","RA","HDL","LDL","Triglycerides","HbA1c")
id_list <- c("ieu-b-5102","ebi-a-GCST90475667","ebi-a-GCST005194","ebi-a-GCST006906",
             "ebi-a-GCST90027158","ebi-a-GCST007090","ebi-a-GCST90018808","ieu-b-7",
             "ebi-a-GCST90038683","ieu-b-41","ebi-a-GCST90014022","ebi-a-GCST90018910",
             "ieu-b-109","ieu-b-110","ieu-b-111","ebi-a-GCST90014006")
continuous_outcomes <- c("BMD", "HDL", "HbA1c", "Triglycerides", "LDL")             
outcome_names <- c(
  `SCZ` = "Schizophrenia",
  `T2D` = "Type 2 Diabetes",
  `stroke` = "stroke",
  `Alzheimer` = "Alzheimer's Disease",
  `Knee` = "Knee Osteoarthritis",
  `Parkinson` = "Parkinson's Disease",
  `Colorectal` = "Colorectal Cancer",
  `IBD` = "Inflammatory Bowel Disease",
  `BD` = "Bipolar Disorder",
  `BMD` = "Bone Mineral Density",
  `RA` = "Rheumatoid Arthritis",
  `CAD` = "Coronary Artery Disease",
  `HDL` = "HDL Cholesterol",
  `LDL` = "LDL Cholesterol",
  `Triglycerides` = "Triglycerides",
  `HbA1c` = "Glycated Hemoglobin")
for (i in 1:length(label)) {
  res_0.05 <- read.csv(paste0(id_list[i],"/results/",label[i],"_MVMR_GRAPPLE_FDR_p_0.05_summary.csv"))
  res_0.1 <- read.csv(paste0(id_list[i],"/results/",label[i],"_MVMR_GRAPPLE_FDR_p_0.1_summary.csv"))
  res_MR <- read.csv(paste0(id_list[i],"/results/",label[i],"_MR_GRAPPLE_FDR_p_0.05_summary.csv"))
  res_0.05$outcome <- label[i]
  res_0.1$outcome <- label[i]
  res_MR$outcome <- label[i]
  all_res_0.05 <- rbind(all_res_0.05,res_0.05)
  all_res_0.1 <- rbind(all_res_0.1,res_0.1)
  all_res_MR <- rbind(all_res_MR,res_MR)
}
all_res_0.05[all_res_0.05$type == "Stepwise","type"] <- "MVMR"
all_res_0.1[all_res_0.1$type == "Stepwise","type"] <- "MVMR"
all_res_MR[all_res_MR$type == "Stepwise","type"] <- "MVMR"
# Use 5e-08, usual ESMR, usual IVW, MRBEE with pleiotropy 0.05
p1 <- all_res_0.05 %>%
  mutate(xintercept = if_else(outcome %in% continuous_outcomes, 0, 1),
                              sig = if_else(pvalue < 0.05/16, "*", "")) %>%
  filter(exposure=="ebi-a-GCST90029070") %>%
  filter(type == "UVMR" | type == "MVMR") %>%
  filter(method != "MRBEE_5e-08_pleio_0") %>%
  filter(method != "GRAPPLE_1e-05") %>%
  filter(method != "IVW_T_5e-08") %>%
  filter(method != "ESMR_optimize_5e-08") %>%
  filter(method != "ESMR_5e-08") %>%
  mutate(method = case_when(method == "IVW_5e-08" ~ "MV-IVW",
                            #method == "ESMR_5e-08" ~ "ESMR",
                            method == "MRBEE_5e-08_pleio_0.05" ~ "MRBEE",
                            method == "GRAPPLE_5e-08" ~ "GRAPPLE",
                            TRUE ~ method)) %>%
  #filter(converge == TRUE | is.na(converge)) %>%
  #filter(se < 1) %>%
  ggplot() +
  geom_vline(aes(xintercept = xintercept)) +
  geom_point(aes(y = type,
                 x = if_else(outcome %in% continuous_outcomes, b, odds),
                 color = method,
                 group = method),
             position=position_dodge(width = 0.9), size = 3) +
  geom_errorbar(aes(y = type,
                    xmin = if_else(outcome %in% continuous_outcomes, log(CI_lower), CI_lower),
                    xmax = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                    color = method),
                position=position_dodge(width = 0.9)) +
  geom_text(aes(y = type,
                x = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                label = sig,
                group = method),
            position = position_dodge(width = 0.9), vjust = -0.5, color = "black",size = 6) +
  xlab("Odds Ratio or Beta Hat (95% CI)") + coord_flip() +
  facet_wrap(~factor(outcome,levels = c('stroke','Knee','CAD','T2D','IBD',
                                        'SCZ','RA','Alzheimer','Parkinson','Colorectal',
                                        'BD','BMD','HDL','LDL','Triglycerides','HbA1c')), ncol=4, scales = "free_y",
             labeller = as_labeller(outcome_names))+
  theme_bw() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=15),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.position = "bottom")
#ggsave("all_summary_0.05.png", plot = p1,width = 15, height = 15)
p2 <- all_res_0.1 %>% mutate(xintercept = if_else(outcome %in% continuous_outcomes, 0, 1)) %>%
  filter(exposure=="ebi-a-GCST90029070") %>%
  filter(type == "UVMR" | type == "MVMR") %>%
  filter(method != "MRBEE_5e-08_pleio_0") %>%
  filter(method != "GRAPPLE_1e-05") %>%
  filter(method != "IVW_T_5e-08") %>%
  filter(method != "ESMR_optimize_5e-08") %>%
  filter(method != "ESMR_5e-08") %>%
  mutate(method = case_when(method == "IVW_5e-08" ~ "MV-IVW",
                            #method == "ESMR_5e-08" ~ "ESMR",
                            method == "MRBEE_5e-08_pleio_0.05" ~ "MRBEE",
                            method == "GRAPPLE_5e-08" ~ "GRAPPLE",
                            TRUE ~ method)) %>%
  #filter(converge == TRUE | is.na(converge)) %>%
  #filter(se < 1) %>%
  ggplot() +
  geom_vline(aes(xintercept = xintercept)) +
  geom_point(aes(y = type,
                 x = if_else(outcome %in% continuous_outcomes, b, odds),
                 color = method,
                 group = method),
             position=position_dodge(width = 0.9), size = 3) +
  geom_errorbar(aes(y = type,
                    xmin = if_else(outcome %in% continuous_outcomes, log(CI_lower), CI_lower),
                    xmax = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                    color = method),
                position=position_dodge(width = 0.9)) +
  xlab("Odds Ratio or Beta Hat (95% CI)") + coord_flip() +
  facet_wrap(~factor(outcome,levels = c('stroke','Knee','CAD','T2D','IBD',
                                        'SCZ','RA','Alzheimer','Parkinson','Colorectal',
                                        'BD','BMD','HDL','LDL','Triglycerides','HbA1c')), ncol=4, scales = "free_y",
             labeller = as_labeller(outcome_names))+
  theme_bw() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=15),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.position = "bottom")
#ggsave("all_summary_0.1.png", plot = p2,width = 15, height = 15)
p3 <- all_res_MR %>% mutate(xintercept = if_else(outcome %in% continuous_outcomes, 0, 1)) %>%
  filter(exposure=="ebi-a-GCST90029070") %>%
  filter(type == "UVMR" | type == "MVMR") %>%
  filter(method != "MRBEE_5e-08_pleio_0") %>%
  filter(method != "GRAPPLE_1e-05") %>%
  filter(method != "IVW_T_5e-08") %>%
  filter(method != "ESMR_optimize_5e-08") %>%
  filter(method != "ESMR_5e-08") %>%
  mutate(method = case_when(method == "IVW_5e-08" ~ "MV-IVW",
                            #method == "ESMR_5e-08" ~ "ESMR",
                            method == "MRBEE_5e-08_pleio_0.05" ~ "MRBEE",
                            method == "GRAPPLE_5e-08" ~ "GRAPPLE",
                            TRUE ~ method)) %>%
  #filter(converge == TRUE | is.na(converge)) %>%
  #filter(se < 1) %>%
  ggplot() +
  geom_vline(aes(xintercept = xintercept)) +
  geom_point(aes(y = type,
                 x = if_else(outcome %in% continuous_outcomes, b, odds),
                 color = method,
                 group = method),
             position=position_dodge(width = 0.9), size = 3) +
  geom_errorbar(aes(y = type,
                    xmin = if_else(outcome %in% continuous_outcomes, log(CI_lower), CI_lower),
                    xmax = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                    color = method),
                position=position_dodge(width = 0.9)) +
  xlab("Odds Ratio or Beta Hat (95% CI)") + coord_flip() +
  facet_wrap(~factor(outcome,levels = c('stroke','Knee','CAD','T2D','IBD',
                                        'SCZ','RA','Alzheimer','Parkinson','Colorectal',
                                        'BD','BMD','HDL','LDL','Triglycerides','HbA1c')), ncol=4, scales = "free_y",
             labeller = as_labeller(outcome_names))+
  theme_bw() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=15),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.position = "bottom")
#ggsave("all_summary_MR.png", plot = p3,width = 15, height = 15)
p4 <- all_res_0.05 %>% mutate(xintercept = if_else(outcome %in% continuous_outcomes, 0, 1)) %>%
  filter(exposure=="ebi-a-GCST90029070") %>%
  filter(type == "only_BMI" | type == "MVMR") %>%
  filter(method != "MRBEE_5e-08_pleio_0") %>%
  filter(method != "GRAPPLE_1e-05") %>%
  filter(method != "IVW_T_5e-08") %>%
  filter(method != "ESMR_optimize_5e-08") %>%
  filter(method != "ESMR_5e-08") %>%
  mutate(method = case_when(method == "IVW_5e-08" ~ "MV-IVW",
                            #method == "ESMR_5e-08" ~ "ESMR",
                            method == "MRBEE_5e-08_pleio_0.05" ~ "MRBEE",
                            method == "GRAPPLE_5e-08" ~ "GRAPPLE",
                            TRUE ~ method)) %>%
  #filter(converge == TRUE | is.na(converge)) %>%
  #filter(se < 1) %>%
  ggplot() +
  geom_vline(aes(xintercept = xintercept)) +
  geom_point(aes(y = type,
                 x = if_else(outcome %in% continuous_outcomes, b, odds),
                 color = method,
                 group = method),
             position=position_dodge(width = 0.9), size = 3) +
  geom_errorbar(aes(y = type,
                    xmin = if_else(outcome %in% continuous_outcomes, log(CI_lower), CI_lower),
                    xmax = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                    color = method),
                position=position_dodge(width = 0.9)) +
  xlab("Odds Ratio or Beta Hat (95% CI)") + coord_flip() +
  facet_wrap(~factor(outcome,levels = c('stroke','Knee','CAD','T2D','IBD',
                                        'SCZ','RA','Alzheimer','Parkinson','Colorectal',
                                        'BD','BMD','HDL','LDL','Triglycerides','HbA1c')), ncol=4, scales = "free_y",
             labeller = as_labeller(outcome_names))+
  theme_bw() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=15),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.position = "bottom")
#ggsave("all_summary_onlyBMI.png", plot = p4,width = 15, height = 15)
p5 <- all_res_0.05 %>%
  mutate(xintercept = if_else(outcome %in% continuous_outcomes, 0, 1),
         sig = case_when(
           pvalue < 0.05/16 ~ "**",
           pvalue < 0.05 ~ "*",
           TRUE ~ ""
         ),
         type = case_when(
           type == "only_BMI" ~ "only BMI",
           TRUE ~ type
         ),
         type = factor(type, levels = c("UVMR", "only BMI", "MVMR"))
  ) %>%
  filter(exposure=="ebi-a-GCST90029070") %>%
  filter(type == "MVMR" | type == "only BMI" | type == "UVMR") %>%
  filter(method == "GRAPPLE_5e-08") %>%
  mutate(method = case_when(method == "IVW_5e-08" ~ "MV-IVW",
                            #method == "ESMR_5e-08" ~ "ESMR",
                            method == "MRBEE_5e-08_pleio_0.05" ~ "MRBEE",
                            method == "GRAPPLE_5e-08" ~ "GRAPPLE",
                            TRUE ~ method)) %>%
  #filter(converge == TRUE | is.na(converge)) %>%
  #filter(se < 1) %>%
  ggplot() +
  geom_vline(aes(xintercept = xintercept)) +
  geom_point(aes(y = type,
                 x = if_else(outcome %in% continuous_outcomes, b, odds),
                 color = type,
                 group = type),
             position=position_dodge(width = 0.9), size = 3) +
  geom_errorbar(aes(y = type,
                    xmin = if_else(outcome %in% continuous_outcomes, log(CI_lower), CI_lower),
                    xmax = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                    color = type),
                position=position_dodge(width = 0.9)) +
  geom_text(aes(y = type,
                x = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                label = sig,
                group = type),
            position = position_dodge(width = 0.9), vjust = -0.05, color = "black",size = 6) +
  xlab("Odds Ratio or Beta Hat (95% CI)") + coord_flip() +
  facet_wrap(~factor(outcome,levels = c('stroke','CAD','HDL','LDL','Triglycerides',
                                        'T2D','HbA1c','IBD','RA','SCZ','BD',
                                        'Alzheimer','Parkinson','Colorectal',
                                        'Knee','BMD')), ncol=4, scales = "free_y",
             labeller = as_labeller(outcome_names))+
  theme_bw() +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 20),
        plot.title = element_text(size=15),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))
ggsave("all_summary_GRAPPLE_New.jpg", plot = p5, device = "jpeg", width = 20, height = 20)
p6 <- all_res_0.05 %>%
  mutate(xintercept = if_else(outcome %in% continuous_outcomes, 0, 1),
         sig = case_when(
           pvalue < 0.05/16 ~ "**",
           pvalue < 0.05 ~ "*",
           TRUE ~ ""
         ),
         type = case_when(
           type == "only_BMI" ~ "only BMI",
           TRUE ~ type
         ),
         type = factor(type, levels = c("UVMR", "only BMI", "MVMR"))
  ) %>%
  filter(exposure=="ebi-a-GCST90029070") %>%
  filter(type == "MVMR" | type == "only BMI" | type == "UVMR") %>%
  filter(method == "IVW_5e-08") %>%
  mutate(method = case_when(method == "IVW_5e-08" ~ "MV-IVW",
                            #method == "ESMR_5e-08" ~ "ESMR",
                            method == "MRBEE_5e-08_pleio_0.05" ~ "MRBEE",
                            method == "GRAPPLE_5e-08" ~ "GRAPPLE",
                            TRUE ~ method)) %>%
  #filter(converge == TRUE | is.na(converge)) %>%
  #filter(se < 1) %>%
  ggplot() +
  geom_vline(aes(xintercept = xintercept)) +
  geom_point(aes(y = type,
                 x = if_else(outcome %in% continuous_outcomes, b, odds),
                 color = type,
                 group = type),
             position=position_dodge(width = 0.9), size = 3) +
  geom_errorbar(aes(y = type,
                    xmin = if_else(outcome %in% continuous_outcomes, log(CI_lower), CI_lower),
                    xmax = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                    color = type),
                position=position_dodge(width = 0.9)) +
  geom_text(aes(y = type,
                x = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                label = sig,
                group = type),
            position = position_dodge(width = 0.9), vjust = -0.05, color = "black",size = 6) +
  xlab("Odds Ratio or Beta Hat (95% CI)") + coord_flip() +
  facet_wrap(~factor(outcome,levels = c('stroke','CAD','HDL','LDL','Triglycerides',
                                        'T2D','HbA1c','IBD','RA','SCZ','BD',
                                        'Alzheimer','Parkinson','Colorectal',
                                        'Knee','BMD')), ncol=4, scales = "free_y",
             labeller = as_labeller(outcome_names))+
  theme_bw() +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 20),
        plot.title = element_text(size=15),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))
ggsave("all_summary_IVW_New.jpg", plot = p6, device = "jpeg", width = 20, height = 20)
p7 <- all_res_0.05 %>%
  mutate(xintercept = if_else(outcome %in% continuous_outcomes, 0, 1),
         sig = case_when(
           pvalue < 0.05/16 ~ "***",
           pvalue < 0.01 ~ "**",
           pvalue < 0.05 ~ "*",
           TRUE ~ ""
         ),
         type = case_when(
           type == "only_BMI" ~ "only BMI",
           TRUE ~ type
         ),
         type = factor(type, levels = c("UVMR", "only BMI", "MVMR"))
  ) %>%
  filter(exposure=="ebi-a-GCST90029070") %>%
  filter(type == "MVMR" | type == "only BMI" | type == "UVMR") %>%
  filter(method == "MRBEE_5e-08_pleio_0.05") %>%
  mutate(method = case_when(method == "IVW_5e-08" ~ "MV-IVW",
                            #method == "ESMR_5e-08" ~ "ESMR",
                            method == "MRBEE_5e-08_pleio_0.05" ~ "MRBEE",
                            method == "GRAPPLE_5e-08" ~ "GRAPPLE",
                            TRUE ~ method)) %>%
  #filter(converge == TRUE | is.na(converge)) %>%
  #filter(se < 1) %>%
  ggplot() +
  geom_vline(aes(xintercept = xintercept)) +
  geom_point(aes(y = type,
                 x = if_else(outcome %in% continuous_outcomes, b, odds),
                 color = type,
                 group = type),
             position=position_dodge(width = 0.9), size = 3) +
  geom_errorbar(aes(y = type,
                    xmin = if_else(outcome %in% continuous_outcomes, log(CI_lower), CI_lower),
                    xmax = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                    color = type),
                position=position_dodge(width = 0.9)) +
  geom_text(aes(y = type,
                x = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                label = sig,
                group = type),
            position = position_dodge(width = 0.9), vjust = -0.05, color = "black",size = 6) +
  xlab("Odds Ratio or Beta Hat (95% CI)") + coord_flip() +
  facet_wrap(~factor(outcome,levels = c('stroke','CAD','HDL','LDL','Triglycerides',
                                        'T2D','HbA1c','IBD','RA','SCZ','BD',
                                        'Alzheimer','Parkinson','Colorectal',
                                        'Knee','BMD')), ncol=4, scales = "free_y",
             labeller = as_labeller(outcome_names))+
  theme_bw() +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 20),
        plot.title = element_text(size=15),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))
ggsave("all_summary_MRBEE_New.jpg", plot = p7, device = "jpeg", width = 20, height = 20)
p8 <- all_res_0.05 %>%
  mutate(xintercept = if_else(outcome %in% continuous_outcomes, 0, 1),
         sig = case_when(
           pvalue < 0.05/16 ~ "***",
           pvalue < 0.01 ~ "**",
           pvalue < 0.05 ~ "*",
           TRUE ~ ""
         ),
         type = case_when(
           type == "only_BMI" ~ "only BMI",
           TRUE ~ type
         ),
         type = factor(type, levels = c("UVMR", "only BMI", "MVMR"))
  ) %>%
  filter(exposure=="ebi-a-GCST90029070") %>%
  filter(type == "MVMR" | type == "only BMI" | type == "UVMR") %>%
  filter(method == "MRBEE_5e-08_pleio_0") %>%
  mutate(method = case_when(method == "IVW_5e-08" ~ "MV-IVW",
                            #method == "ESMR_5e-08" ~ "ESMR",
                            method == "MRBEE_5e-08_pleio_0" ~ "MRBEE",
                            method == "GRAPPLE_5e-08" ~ "GRAPPLE",
                            TRUE ~ method)) %>%
  #filter(converge == TRUE | is.na(converge)) %>%
  #filter(se < 1) %>%
  ggplot() +
  geom_vline(aes(xintercept = xintercept)) +
  geom_point(aes(y = type,
                 x = if_else(outcome %in% continuous_outcomes, b, odds),
                 color = type,
                 group = type),
             position=position_dodge(width = 0.9), size = 3) +
  geom_errorbar(aes(y = type,
                    xmin = if_else(outcome %in% continuous_outcomes, log(CI_lower), CI_lower),
                    xmax = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                    color = type),
                position=position_dodge(width = 0.9)) +
  geom_text(aes(y = type,
                x = if_else(outcome %in% continuous_outcomes, log(CI_higher), CI_higher),
                label = sig,
                group = type),
            position = position_dodge(width = 0.9), vjust = -0.05, color = "black",size = 6) +
  xlab("Odds Ratio or Beta Hat (95% CI)") + coord_flip() +
  facet_wrap(~factor(outcome,levels = c('stroke','CAD','HDL','LDL','Triglycerides',
                                        'T2D','HbA1c','IBD','RA','SCZ','BD',
                                        'Alzheimer','Parkinson','Colorectal',
                                        'Knee','BMD')), ncol=4, scales = "free_y",
             labeller = as_labeller(outcome_names))+
  theme_bw() +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 20),
        plot.title = element_text(size=15),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))
ggsave("all_summary_MRBEE_New_nopleio.jpg", plot = p8, device = "jpeg", width = 20, height = 20)

