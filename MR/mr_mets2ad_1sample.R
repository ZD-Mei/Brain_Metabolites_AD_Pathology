rm(list = ls())
library(tidyverse)
library(data.table)
library(openxlsx)
library(sjmisc)
library(TwoSampleMR)
library(biomaRt)
library(ieugwasr)


setwd("workpath")

select <- dplyr::select

load("data/dat_combine_all.RData")

mets_list <- read_lines("lists/ad_path_mets_list.txt")

# Extract SNPs ------
met <- "METID"
(met_name <- brain_targeted_info$Compounds[brain_targeted_info$com_id==met])

re_clump_filtered <- 
  read.xlsx(paste0("gwas/brainMetabGWAS/resClump/gwas_",
                   met,"_",met_name,"_clump_filtered.xlsx"))
(snps <- paste0("SNP.",re_clump_filtered$chr.exposure,
                "_",re_clump_filtered$pos.exposure))
(snpIDs <- paste(re_clump_filtered$SNP,collapse = ";"))
(nsnps <- length(re_clump_filtered$SNP))


## IV data -----
dat_iv <- 
  fread("gwas/data/mets_sig_snps.raw")
head(dat_iv);dim(dat_iv) 
dat_iv <- dat_iv %>% 
  select(-c(FID,PAT,MAT,SEX,PHENOTYPE))
colnames(dat_iv) <- gsub(":","_",colnames(dat_iv))
colnames(dat_iv)[2:ncol(dat_iv)] <- 
  paste0("SNP.",colnames(dat_iv)[2:ncol(dat_iv)])
colnames(dat_iv) <- gsub("_[^_]+$","",colnames(dat_iv))

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
  select(FID,IID,age_death,msex,pmi,amylsqrt,tangsqrt,cogng_demog_slope,
         brain1_com1:brain1_com177)
table(colSums(is.na(dat_pheno)))
head(dat_pheno);dim(dat_pheno)

## data combine -----
dat_com <- dat_pc2 %>% 
  left_join(dat_iv,by="IID") %>% 
  inner_join(dat_pheno,by=c("FID","IID"))
head(dat_com);dim(dat_com) 
colSums(is.na(dat_com))

# onesample MR amyloid ------
library(AER)
table(snps %in% colnames(dat_com))

covs <- "age_death + msex + pmi + PC1 + PC2 + PC3 + PC4 + BU + Broad"

## amyloid 
iv_form_amy <- 
  as.formula(paste0("amylsqrt ~ ",met,"+",covs,"|",
                    paste(snps,collapse = "+"),"+",covs))
iv_mod_amy <- ivreg(iv_form_amy,data = dat_com)

sum_iv_amy <- summary(iv_mod_amy, diagnostics = TRUE)
F_stat_amy <- sum_iv_amy$diagnostics[grep("Weak", rownames(sum_iv_amy$diagnostics)),
                             "statistic"]
F_stat_p_amy <- sum_iv_amy$diagnostics[grep("Weak", rownames(sum_iv_amy$diagnostics)),
                               "p-value"]
re_amy_1mr <- c(met,"amyloid",
                coef(summary(iv_mod_amy))[2,c(1,2,4)],
                nsnps,snpIDs,F_stat_amy,F_stat_p_amy)

## tangles 
iv_form_tang <- 
  as.formula(paste0("tangsqrt ~ ",met,"+",covs,"|",
                    paste(snps,collapse = "+"),"+",covs))
iv_mod_tang <- ivreg(iv_form_tang,data = dat_com)
sum_iv_tang <- summary(iv_mod_tang, diagnostics = TRUE)
F_stat_tang <- sum_iv_tang$diagnostics[grep("Weak", rownames(sum_iv_tang$diagnostics)),
                                     "statistic"]
F_stat_p_tang <- sum_iv_tang$diagnostics[grep("Weak", rownames(sum_iv_tang$diagnostics)),
                                       "p-value"]
re_tang_1mr <- c(met,"tangles",
                 coef(summary(iv_mod_tang))[2,c(1,2,4)],
                 nsnps,snpIDs,F_stat_tang,F_stat_p_tang)

## cognitive slope
iv_form_cog <- 
  as.formula(paste0("cogng_demog_slope ~ ",met,"+",covs,"|",
                    paste(snps,collapse = "+"),"+",covs))
iv_mod_cog <- ivreg(iv_form_cog,data = dat_com)
sum_iv_cog <- summary(iv_mod_cog, diagnostics = TRUE)
F_stat_cog <- sum_iv_cog$diagnostics[grep("Weak", rownames(sum_iv_cog$diagnostics)),
                                     "statistic"]
F_stat_p_cog <- sum_iv_cog$diagnostics[grep("Weak", rownames(sum_iv_cog$diagnostics)),
                                       "p-value"]
re_cog_1mr <- c(met,"cognitive",
                coef(summary(iv_mod_cog))[2,c(1,2,4)],
                nsnps,snpIDs,F_stat_cog,F_stat_p_cog)


## rbind results 
re_1mr <- rbind(re_amy_1mr,re_tang_1mr,re_cog_1mr) %>% 
  as.data.frame(.)
re_1mr
colnames(re_1mr) <- 
  c("com_id","outcome","Beta","SE","P_val",
    "nsnps","snpIDs","F_statistic","F_statistic_pval")
re_1mr[,c(3:6,8:9)] <- map_df(re_1mr[,c(3:6,8:9)],as.numeric)
fwrite(re_1mr,paste0("results/mr/mets2ad_1sample/",met,"_to_ad_1sample.csv"))



