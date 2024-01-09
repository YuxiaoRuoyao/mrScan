library(TwoSampleMR)
library(dplyr)
library(data.table)
library(GRAPPLE)

ex_dat1 <- readRDS(snakemake@input[["file1"]]) # X/Y + M
ex_dat2 <- readRDS(snakemake@input[["file2"]]) # Z
ex_dat3 <- readRDS(snakemake@input[["file3"]]) # Z + M
ex_dat4 <- readRDS(snakemake@input[["file4"]]) # X/Y
out <- snakemake@output[["out"]]


ID1 <- unique(ex_dat4$id.exposure) # X/Y
ID2 <- unique(ex_dat2$id.exposure) # Z
ID3 <- unique(ex_dat1$id.exposure)[-1] # M
# Calculate cor between X/Y and Z
out_dat1 <- extract_outcome_data(snps = ex_dat4$SNP,outcomes = ID2)
out_dat2 <- extract_outcome_data(snps = ex_dat2$SNP,outcomes = ID1)
dat_1_2 <- harmonise_data(ex_dat1, out_dat1)
dat_2_1 <- harmonise_data(ex_dat2, out_dat2)
X_1_2 <- dat_1_2 %>%
  rename(beta1 = beta.exposure, beta2 = beta.outcome) %>%
  select(SNP, beta1, beta2)
X_2_1 <- dat_2_1 %>%
  rename(beta2 = beta.exposure, beta1 = beta.outcome) %>%
  select(SNP, beta1, beta2) %>% filter(!SNP %in% X_1_2$SNP)
X <- bind_rows(X_1_2, X_2_1)
cor_vals<- data.frame(id1 = ID1,id2 = ID2, cor = with(X, cor(beta1, beta2)))
# MVMR by GRAPPLE
out_1_2 <- extract_outcome_data(snps = ex_dat1$SNP,outcomes = ID2)
out_3_4 <- extract_outcome_data(snps = ex_dat3$SNP,outcomes = ID1)
mvdat_1 <- mv_harmonise_data(ex_dat1,out_1_2)
mvdat_2 <- mv_harmonise_data(ex_dat3,out_3_4)

grapple.dat1<-cbind(mvdat_1$exposure_beta,mvdat_1$exposure_se,
                    mvdat_1$outcome_beta,mvdat_1$outcome_se)
grapple.dat2<-cbind(mvdat_2$exposure_beta,mvdat_2$exposure_se,
                    mvdat_2$outcome_beta,mvdat_2$outcome_se)
i <- ncol(mvdat_1$exposure_beta)
j <- ncol(mvdat_2$exposure_beta)
colnames(grapple.dat1)<-c(paste0("gamma_exp",seq(1,i)),
                          paste0("se_exp",seq(1,i)),"gamma_out1","se_out1")
colnames(grapple.dat2)<-c(paste0("gamma_exp",seq(1,j)),
                          paste0("se_exp",seq(1,j)),"gamma_out1","se_out1")
grapple.res1<-grappleRobustEst(grapple.dat1,plot.it =FALSE,niter = 100000)
grapple.res2<-grappleRobustEst(grapple.dat2,plot.it =FALSE,niter = 100000)
res1 <- data.frame(id.exposure=colnames(mvdat_1$exposure_beta),
                   id.outcome=mvdat_1$outname$id.outcome,
                   b=grapple.res1$beta.hat,
                   se=sqrt(diag(grapple.res1$beta.var)),
                   pval=grapple.res1$beta.p.value) %>%
        filter(id.exposure == ID1)
res2 <- data.frame(id.exposure=colnames(mvdat_2$exposure_beta),
                   id.outcome=mvdat_2$outname$id.outcome,
                   b=grapple.res2$beta.hat,
                   se=sqrt(diag(grapple.res2$beta.var)),
                   pval=grapple.res2$beta.p.value) %>%
        filter(id.exposure == ID2)
saveRDS(list(mr12 = res1, mr21 = res2, cor = cor_vals),
        file = out)
