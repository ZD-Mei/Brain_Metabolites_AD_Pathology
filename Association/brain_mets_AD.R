rm(list = ls())
library(tidyverse)
library(data.table)
library(openxlsx)
library(sjmisc)
library(grid)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(ComplexUpset)
library(RColorBrewer)
library(ggsci)
library(scales)
library(ggrepel)


setwd("workpath")

select <- dplyr::select

load("data/dat_combine_all.RData")

dat_brain_nona <- dat_brain %>% 
  filter(!is.na(brain1_com1))

# association -------
outcome_con <- c("amylsqrt","tangsqrt","cogng_demog_slope") 
adj_vars <- c("cogng_demog_slope")

(i=brain_targeted_list[2])
(j=outcome_con[2])

re_brain_con <- c()
for (i in brain_targeted_list) {
  for (j in outcome_con) {
    if (i %in% adj_vars){
      mod.form <- as.formula(paste0(j,"~",i,"+age_death+pmi"))
    } else {
      mod.form <- as.formula(paste0(j,"~",i,"+age_death+msex+educ+pmi")) 
    }
    mod.log <- lm(mod.form,data = dat_all)
    re_brain_con <- rbind(re_brain_con,c(i,j,coef(summary(mod.log))[2,c(1,2,4)]))
  }
}
re_brain_con <- as.data.frame(re_brain_con)
head(re_brain_con);dim(re_brain_con) #531
colnames(re_brain_con) <- c("com_id","outcome","Beta","SE","P")
re_brain_con[,3:5] <- map_df(re_brain_con[,3:5],as.numeric)

re_brain_con2 <- re_brain_con %>% 
  inner_join(brain_targeted_info,by="com_id") %>% 
  group_by(outcome) %>% 
  mutate(fdr=p.adjust(P, n = length(P),method = "BH")) %>% 
  ungroup(.)
head(re_brain_con2) 
table(re_brain_con2$P < 0.05)  
table(re_brain_con2$fdr < 0.1) 
# write.xlsx(re_brain_con2,"results/brain_mets_outcome.xlsx")

re_amyloid_sig <- re_brain_con2 %>% 
  filter(outcome=="amylsqrt" & fdr<0.05)
dim(re_amyloid_sig)

re_amyloid_sig <- re_brain_con2 %>% 
  filter(outcome=="tangsqrt" & fdr<0.05)
dim(re_amyloid_sig)

# output results -----
head(re_brain_con2);dim(re_brain_con2)
write.xlsx(re_brain_con2,"results/brain_mets_outcome/brain_mets_outcome.xlsx")


# output for STable -----
table(re_brain_con2$outcome)

re_brain_ad_cog_wide <- re_brain_con2 %>% 
  select(-c("HMDB_main")) %>% #"com_id",
  pivot_wider(values_from = c(Beta,SE,P,fdr),
              names_from = outcome,
              names_vary = "slowest") %>% 
  select(com_id,Compounds,Class,everything())
head(re_brain_ad_cog_wide)
write.xlsx(re_brain_ad_cog_wide,
           "results/brain_mets_outcome/stable_brain_mets_ad_cog.xlsx")

