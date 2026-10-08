rm(list = ls())

library(data.table)
library(tidyverse)
library(glmnet)
library(gtools)
library(openxlsx)
library(broom)
library(table1)

setwd("workpath")

select <- dplyr::select
int <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}

# GWAS data ------
gwas_path <- "Imputation_HRC_v2.3/"

dat_pc <- fread(paste0(gwas_path,"meta_data.txt"),
                colClasses = c(FID = "character", IID = "character"))
head(dat_pc);dim(dat_pc) 
dat_pc$projid <- as.numeric(dat_pc$FID)

# snRNAseq data ------
snRNAseq <- 
  read_rds("pastel_resources_data2share_snRNA_pseudoBulk_7majCellTypes.rds")
head(snRNAseq$ext[,1:10]);dim(snRNAseq$ext) 

ext <- snRNAseq$ext %>% 
  rotate_df(.) %>% 
  rownames_to_column(var = "projid") %>% 
  mutate(projid=as.numeric(projid))
head(ext[,1:5]);dim(ext) #424/16835

# metabolites data -----
load("data/dat_combine_all.RData")
head(colnames(dat_brain),80);dim(dat_brain) # 
table(is.na(dat_brain$brain1_com1))

dat_table1 <- dat_brain[,1:77]
colnames(dat_table1)

table(dat_table1$projid %in% dat_pc$projid)

dat_table1 <- dat_table1 %>% 
  mutate(GWAS=ifelse(projid %in% dat_pc$projid,"Participants with genotyping data","No")) %>% 
  mutate(snRNA=ifelse(projid %in% ext$projid,"Participants with sn-RNA seq data","No"))
table(dat_table1$snRNA)
table(dat_table1$GWAS)
summary(dat_table1$age_death - dat_table1$age_at_visit)

# clean variables -------
## factors ----- 
table(dat_table1$msex)
dat_table1$msex <- 
  factor(dat_table1$msex, levels=c(0,1),
         labels=c("Female","Male"))

table(dat_table1$race7)
dat_table1$race7 <- 
  factor(dat_table1$race7, levels=c(1,2,3),
         labels=c("White","Black","American Indian/Alaska Native"))

table(dat_table1$apoe_genotype)
table(dat_table1$apoe4_b2)
dat_table1$apoe4_b2 <- 
  factor(dat_table1$apoe4_b2,levels = c(0,1),
         labels = c("Non-carrier","Carrier"))


## set names -----
label(dat_table1$age_bl) <- "Age at baseline"
label(dat_table1$age_death) <- "Age at death"
label(dat_table1$age_at_visit) <- "Age at blood sample collection"
label(dat_table1$apoe4_b2) <- "APOE ε4"
label(dat_table1$msex) <- "Sex"
label(dat_table1$race7) <- "Race"
label(dat_table1$educ) <- "Education"
label(dat_table1$amylsqrt) <- "Amyloid-β burden"
label(dat_table1$tangsqrt) <- "Tau tangle density"
label(dat_table1$cogng_demog_slope) <- "Cognitive slope"
label(dat_table1$pmi) <- "Post-mortem interval"
label(dat_table1$bmi) <- "BMI at blood sample collection"
label(dat_table1$phys5itemsum) <- "Physical activity at blood sample collection"

## set units -----
units(dat_table1$age_bl) <- "years"
units(dat_table1$age_death) <- "years"
units(dat_table1$age_at_visit) <- "years"
units(dat_table1$educ) <- "years"
units(dat_table1$pmi) <- "hours"
units(dat_table1$bmi) <- "kg/m2"
units(dat_table1$phys5itemsum) <- "hours/week"

# table1  -----
table1(~age_bl+msex+race7+apoe4_b2+educ+bmi+phys5itemsum+age_death+pmi
       +amylsqrt+tangsqrt+cogng_demog_slope|snRNA,
       overall = c(left="Overall"),
       render.continuous=c(.="Mean (SD)"),
       # render.missing = NULL,
       data=dat_table1)

table1(~age_bl+msex+race7+apoe4_b2+educ+bmi+phys5itemsum+age_death+pmi
       +amylsqrt+tangsqrt+cogng_demog_slope|GWAS,
       overall = c(left="Overall"),
       render.continuous=c(.="Mean (SD)"),
       # render.missing = NULL,
       data=dat_table1)

