rm(list = ls())
library(tidyverse)
library(data.table)
library(openxlsx)
library(sjmisc)
library(grid)
library(gridExtra)
library(missRanger)
library(pheatmap)
library(RColorBrewer)
library(ComplexUpset)
library(RColorBrewer)
library(ggsci)
library(scales)
library(ggrepel)
library(nnet)
library(paletteer)
library(TwoSampleMR)
select <- dplyr::select

setwd("workpath")

load("data/dat_combine_all.RData")

# mets list --------
mets_list <- read_lines("lists/ad_path_mets_list.txt")
head(mets_list);length(mets_list)


## IV data -----
dat_iv <- fread("gwas/data/ad_sig_snps.raw")
head(dat_iv);dim(dat_iv)
dat_iv <- dat_iv %>% 
  select(IID,`17:8833591_C`,`19:45411941_T`)
colnames(dat_iv) <- gsub(":","_",colnames(dat_iv))
colnames(dat_iv)[2:3] <- paste0("SNP.",colnames(dat_iv)[2:3])

## PC data -----
dat_pc <- 
  fread("/data/junli-lab/RUSH_Data/GWAS/Imputation_HRC_v2.3/meta_data.txt",
        colClasses = c(FID = "character", IID = "character"))
head(dat_pc);dim(dat_pc) #3383
table(dat_pc$rawdata_source)
dat_pc2 <- dat_pc %>% 
  mutate(FID = 0,
         IID=paste0(IID,"_",IID),
         BU=case_when(
           rawdata_source=="BU" ~ 1,
           TRUE ~ 0),
         Broad=case_when(
           rawdata_source=="BROAD" ~ 1,
           TRUE ~ 0)) %>%
  select(FID,IID,PC1:PC4,BU,Broad)
head(dat_pc2)

## Phenotype data -----
colnames(dat_brain)

dat_pheno <- dat_brain %>% 
  mutate(projid_update=str_pad(projid,width = 8,pad = "0",side = "left"),
         FID=0,
         IID=paste0(projid_update,"_",projid_update)
  ) %>% 
  select(FID,IID,age_death,msex,pmi,amylsqrt,tangsqrt,brain1_com1:brain1_com177)
table(colSums(is.na(dat_pheno)))
head(dat_pheno);dim(dat_pheno)

## data combine -----
dat_com <- dat_pc2 %>% 
  left_join(dat_iv,by="IID") %>% 
  inner_join(dat_pheno,by=c("FID","IID")) #%>% 
head(dat_com);dim(dat_com)
colSums(is.na(dat_com))

## onesample MR tangle ------
library(AER)
covs <- "age_death + msex + pmi + PC1 + PC2 + PC3 + PC4 + BU + Broad"
(mets <- mets_list[1])

re_tang_1mr <- c()
for (mets in mets_list) {
  iv_form <- as.formula(paste0(mets,"~ tangsqrt +",covs,"|",
                               "SNP.19_45411941_T + SNP.17_8833591_C + ",covs))
  iv_mod <- ivreg(iv_form,data = dat_com)
  sum_iv <- summary(iv_mod, diagnostics = TRUE)
  F_stat <- sum_iv$diagnostics[grep("Weak", rownames(sum_iv$diagnostics)),
                               "statistic"]
  F_stat_p <- sum_iv$diagnostics[grep("Weak", rownames(sum_iv$diagnostics)),
                                 "p-value"]
  re_tang_1mr <- 
    rbind(re_tang_1mr,c(mets,"tangle",coef(summary(iv_mod))[2,c(1,2,4)],
                       F_stat,F_stat_p))
}
re_tang_1mr <- as.data.frame(re_tang_1mr)
head(re_tang_1mr);dim(re_tang_1mr) #46
colnames(re_tang_1mr) <- c("com_id","exposure","Beta","SE","P_val",
                           "F_statitic","F_statistic_pval")
re_tang_1mr_info <- re_tang_1mr %>% 
  left_join(brain_targeted_info,by="com_id") %>% 
  arrange(P_val)
head(re_tang_1mr_info)
re_tang_1mr_info[,3:7] <- map_df(re_tang_1mr_info[,3:7],as.numeric)
write.xlsx(re_tang_1mr_info,"results/mr/mr_tangle_to_mets_1sample.xlsx")


dat_com_nona <- dat_com %>% drop_na(.)
stage1 <- lm(tangsqrt ~ SNP.19_45411941_T + SNP.17_8833591_C + age_death + msex + pmi,
             data = dat_com_nona)
summary(stage1)
x_hat <- fitted(stage1)  

stage2 <- lm(dat_com_nona$brain1_com75 ~ x_hat + dat_com_nona$age_death
             + dat_com_nona$msex + dat_com_nona$pmi)
summary(stage2)

iv_mod <- ivreg(brain1_com75 ~ tangsqrt + age_death + msex + pmi | 
                  SNP.19_45411941_T + SNP.17_8833591_C  + age_death + msex + pmi,
                data = dat_com_nona)
summary(iv_mod)
coef(summary(iv_mod))[2,c(1,2,4)]

## onesample MR amyloid ------
covs <- "age_death + msex + pmi + PC1 + PC2 + PC3 + PC4 + BU + Broad"
(mets <- mets_list[1])
re_amy_1mr <- c()
for (mets in mets_list) {
  iv_form <- as.formula(paste0(mets,"~ amylsqrt +",covs,"|",
                               "SNP.19_45411941_T + ",covs))
  iv_mod <- ivreg(iv_form,data = dat_com)
  sum_iv <- summary(iv_mod, diagnostics = TRUE)
  F_stat <- sum_iv$diagnostics[grep("Weak", rownames(sum_iv$diagnostics)),
                               "statistic"]
  F_stat_p <- sum_iv$diagnostics[grep("Weak", rownames(sum_iv$diagnostics)),
                                 "p-value"]
  re_amy_1mr <- 
    rbind(re_amy_1mr,c(mets,"amyloid",coef(summary(iv_mod))[2,c(1,2,4)],
                       F_stat,F_stat_p))
}
re_amy_1mr <- as.data.frame(re_amy_1mr)
head(re_amy_1mr);dim(re_amy_1mr) #63
colnames(re_amy_1mr) <- c("com_id","exposure","Beta","SE","P_val",
                          "F_statitic","F_statistic_pval")
re_amy_1mr_info <- re_amy_1mr %>% 
  left_join(brain_targeted_info,by="com_id") %>% 
  arrange(P_val)
head(re_amy_1mr_info)
re_amy_1mr_info[,3:7] <- map_df(re_amy_1mr_info[,3:7],as.numeric)
write.xlsx(re_amy_1mr_info,"results/mr/mr_amyloid_to_mets_1sample.xlsx")



