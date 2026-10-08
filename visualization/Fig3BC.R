rm(list = ls())
library(tidyverse)
library(data.table)
library(grid)
library(gridExtra)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)
library(sjmisc)
library(ComplexUpset)
library(RColorBrewer)
library(ggsci)
library(scales)
library(ggrepel)


setwd("workpath")
rushpath <- "rushpath"

select <- dplyr::select


load("data/dat_combine_all.RData")


# association results -----
re_assoc_all <- read.xlsx("results/brain_mets_outcome/brain_mets_outcome.xlsx")
table(re_assoc_all$outcome)
re_brain_assoc <- re_assoc_all %>% 
  filter(outcome %in% c("amylsqrt","tangsqrt","cogng_demog_slope"))

re_amy_sig <- re_assoc_all %>% 
  filter(outcome == "amylsqrt" & fdr <0.05)
dim(re_amy_sig)

re_tang_sig <- re_assoc_all %>% 
  filter(outcome == "tangsqrt" & fdr <0.05)
dim(re_tang_sig)


# APOE4 data--------
table(dat_brain$apoe4_b2)
dat_brain_ap4 <- dat_brain %>%
  filter(apoe4_b2==1)
dim(dat_brain_ap4) 

dat_brain_nap4 <- dat_brain %>%
  filter(apoe4_b2==0)
dim(dat_brain_nap4) 


## Continuous outcome --------
run_lm <- function(data, mets_list, outcome_list, adj_list) {
  results <- c()
  for (i in mets_list) {
    for (j in outcome_list) {
      # Adjust the formula based on whether i is in the adj_vars
      if (i %in% adj_list) {
        mod.form <- as.formula(paste0(j, "~", i, "+age_death+pmi"))
      } else {
        mod.form <- as.formula(paste0(j, "~", i, "+age_death+msex+educ+pmi"))
      }
      mod.lm <- lm(mod.form, data = data)
      results <- rbind(results, c(i, j, coef(summary(mod.lm))[2, c(1, 2, 4)]))
    }
  }
  results_df <- as.data.frame(results)
  colnames(results_df) <- c("com_id", "outcome", "Beta", "SE", "P")
  results_df[, 3:5] <- map_df(results_df[, 3:5], as.numeric)
  return(results_df)
}

# Continuous outcomes
outcome_con <- c("amylsqrt","tangsqrt","cogng_demog_slope")

# Adjustment variables
adj_vars <- c("cogng_demog_slope")

### APOE4 carrier ---------
re_brain_con_ap4 <- 
  run_lm(dat_brain_ap4, brain_targeted_list, outcome_con, adj_vars)

# Adjust for FDR and merge with brain_targeted_info
re_brain_con_ap42 <- re_brain_con_ap4 %>%
  inner_join(brain_targeted_info, by = "com_id") %>%
  group_by(outcome) %>%
  mutate(fdr = p.adjust(P, method = "BH")) %>%
  ungroup()

# Check results
head(re_brain_con_ap42)
table(re_brain_con_ap42$P < 0.05) # 133
table(re_brain_con_ap42$fdr < 0.05) # 40

### APOE4 non-carrier ---------
re_brain_con_nap4 <- 
  run_lm(dat_brain_nap4, brain_targeted_list, outcome_con, adj_vars)

# Adjust for FDR and merge with brain_targeted_info
re_brain_con_nap42 <- re_brain_con_nap4 %>%
  inner_join(brain_targeted_info, by = "com_id") %>%
  group_by(outcome) %>%
  mutate(fdr = p.adjust(P, method = "BH")) %>%
  ungroup()

# Check results
head(re_brain_con_nap42)
table(re_brain_con_nap42$P < 0.05) #153
table(re_brain_con_nap42$fdr < 0.05) #99

# _Strata results combine -------
re_brain_total <- read.xlsx("results/brain_mets_outcome/brain_mets_outcome.xlsx")
head(re_brain_total);dim(re_brain_total) #531
re_brain_total$Data = "Total"
table(re_brain_total$outcome)
re_brain_total <- re_brain_total %>% 
  filter(outcome %in% outcome_con)

re_brain_ap4 <- re_brain_con_ap42 %>% 
  mutate(Data="APOE4 carriers")
head(re_brain_ap4);dim(re_brain_ap4) #531/10

re_brain_nap4 <- re_brain_con_nap42 %>% 
  mutate(Data="APOE4 non-carriers")
head(re_brain_nap4);dim(re_brain_nap4) #531/10

re_brain_all <- bind_rows(re_brain_total,re_brain_ap4,re_brain_nap4)
head(re_brain_all);dim(re_brain_all) #1593
re_brain_all$pair <- paste0(re_brain_all$com_id,"_",re_brain_all$outcome)
write.xlsx(re_brain_all,"results/strata_apoe/brain_mets_outcome_withstrata.xlsx")

## Interaction results -----
re_brain_int <- read.xlsx("results/interaction_apoe/brain_mets_apoe4_int_AD.xlsx")
head(re_brain_int);dim(re_brain_int)
re_brain_int <- re_brain_int %>% 
  select(-c(vars,Compounds,Class,HMDB_main))
colnames(re_brain_int)[3:6] <- paste0(colnames(re_brain_int[3:6]),"_int")

re_brain_int_sig <- re_brain_int %>% 
  filter(fdr_int <0.05) %>% 
  mutate(pair1=paste0(com_id,"_",outcome))
length(unique(re_brain_int_sig$com_id)) #22
re_brain_int_sig$pair1

re_brain_int_sig_output <- re_brain_int_sig %>% 
  filter(!duplicated(com_id))
dim(re_brain_int_sig_output)
write.xlsx(re_brain_int_sig_output,
           "results/apoe4_interaction_sig_metabolites.xlsx")

re_int_sig <- re_brain_int_sig %>% 
  group_by(outcome) %>% 
  summarise(sig_n=n()) %>% 
  arrange(desc(sig_n))
re_int_sig

re_int_amy <- re_brain_int %>% 
  filter(outcome=="amylsqrt") %>% 
  filter(com_id %in% re_amy_sig$com_id) %>% 
  mutate(fdr_int_amy=p.adjust(P_int,method = "BH")) %>% 
  arrange(fdr_int_amy)
re_int_amy

re_int_tang <- re_brain_int %>% 
  filter(outcome=="tangsqrt") %>% 
  filter(com_id %in% re_tang_sig$com_id) %>% 
  mutate(fdr_int_tang=p.adjust(P_int,method = "BH")) %>% 
  arrange(fdr_int_tang)
re_int_tang

# output for Stable -----
re_strata <- re_brain_all %>% 
  filter(outcome %in% c("amylsqrt","tangsqrt")) %>% 
  filter(Data != "Total") %>% 
  select(-pair)
head(re_strata)
table(re_strata$outcome)

re_strata_wide1 <- re_strata %>% 
  pivot_wider(values_from = c(Beta,SE,P,fdr),
              names_from = Data,
              names_vary = "slowest")
head(re_strata_wide1)
re_strata_wide2 <- re_strata_wide1 %>% 
  left_join(re_brain_int,by=c("com_id","outcome")) %>% 
  pivot_wider(values_from = `Beta_APOE4 carriers`:fdr_int,
              names_from = outcome,
              names_vary = "slowest") %>% 
  select(-c(com_id,HMDB_main))
head(re_strata_wide2)
write.xlsx(re_strata_wide2,"results/strata_apoe/stable_strata_all.xlsx")

re_strata_amy <- re_strata_wide1 %>%
  filter(outcome=="amylsqrt" & com_id %in% re_amy_sig$com_id) %>% 
  left_join(re_int_amy,by=c("com_id","outcome")) %>% 
  select(-c(com_id,HMDB_main))
head(re_strata_amy);dim(re_strata_amy)
write.xlsx(re_strata_amy,"results/strata_apoe/stable_strata_amyloid.xlsx")

re_strata_tang <- re_strata_wide1 %>%
  filter(outcome=="tangsqrt" & com_id %in% re_tang_sig$com_id) %>% 
  left_join(re_int_tang,by=c("com_id","outcome")) %>% 
  select(-c(com_id,HMDB_main))
head(re_strata_tang);dim(re_strata_tang) 
write.xlsx(re_strata_tang,"results/strata_apoe/stable_strata_tangle.xlsx")


## combine ------
re_tang_comb <- re_brain_all %>% 
  filter(outcome=="tangsqrt",com_id %in% re_tang_sig$com_id) %>% 
  inner_join(re_int_tang,by=c("com_id","outcome")) %>% 
  filter(fdr_int_tang<0.05)
dim(re_tang_comb) #75

re_tang_comb <- re_tang_comb %>% 
  mutate(
    sig_ass=case_when(
      P>=0.05 ~ "",
      P<0.05 & fdr >=0.05 ~ "*",
      fdr <0.05 ~ '#'),
    sig_int=case_when(
      P_int>=0.05 ~ "",
      P_int<0.05 & fdr_int_tang >=0.05 ~ "*",
      fdr_int_tang <0.05 ~ '#')
  ) 
head(re_tang_comb);dim(re_tang_comb) # 75


re_tang_comb$Data <- 
  factor(re_tang_comb$Data,
         levels = c("Total","APOE4 carriers","APOE4 non-carriers"))

re_tang_comb <- re_tang_comb %>%
  filter(Data != "Total") %>% 
  mutate(Beta_max = max(Beta+2.2*SE, na.rm = TRUE),
         Beta_min = min(Beta-2*SE, na.rm = TRUE))
head(re_tang_comb)

## forest plot -------
library(ggsci)
library(paletteer)
paletteer_d("ggsci::default_jama")

(col_pal=pal_nejm("default",alpha=0.8)(8))
scales::show_col(col_pal)
my_color <- c("Acylcarnitines" = col_pal[1],
              "Aminoacids" = col_pal[2],
              "Biogenic Amines" = col_pal[3],
              "Glycerides" = col_pal[4],
              "Glycerophospholipids" = col_pal[5],
              "Sphingolipids" = col_pal[6],
              "Cholesterol Esters" = col_pal[7],
              "Sugars" = col_pal[8])

table(re_tang_comb$Class)
re_tang_comb$Class <- 
  factor(re_tang_comb$Class,
         levels = c("Aminoacids","Biogenic Amines","Acylcarnitines",
                    "Sphingolipids","Glycerophospholipids"))

re_tang_comb <- re_tang_comb %>% 
  arrange(Class,Compounds)

met_levels <- c(
  setdiff(unique(re_tang_comb$Compounds), c("LPC(18:1)", "LPC-O(18:1)")),
  "LPC(18:1)",
  "LPC-O(18:1)"
)

re_tang_comb$Compounds <- 
  factor(re_tang_comb$Compounds,
         levels = met_levels)

cols <- c("#084081", "#7bccc4")
cols <- c("#0F6673","#61B98C")

p_fp <- ggplot(re_tang_comb,aes(x=Compounds,y=Beta,color=Data))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_point(position = position_dodge(0.8),size=3)+
  geom_errorbar(aes(ymin= Beta-1.96*SE, ymax = Beta+1.96*SE),
                width=0,position = position_dodge(0.8), size = 1)+
  scale_color_manual(values = c("#c96353","#338ec4"),name="")+
  labs(x="Compounds",y="Beta (95% CI) for Tau tangles")+
  theme_bw()+
  theme(text=element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
p_fp

p_fp_lgd <- p_fp+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.position = "right")
p_fp_legend <- ggpubr::get_legend(p_fp_lgd) 

pdf("figures/strata_by_apoe/brain_mets_strata_legend.pdf", height = 1,width = 2)
grid.draw(p_fp_legend)
dev.off()

dat_class <- re_tang_comb
p_class <- ggplot(dat_class,aes(x=Compounds,y=1,fill=Class))+
  geom_col(width = 1)+
  coord_cartesian(expand = FALSE)+
  scale_fill_manual(values = my_color)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,color="black"),
        text=element_text(size = 20),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        axis.ticks.y = element_blank())
p_class

p_class_lgd <- p_class+
  theme(legend.position = "right")
p_class_legend <- ggpubr::get_legend(p_class_lgd) 

pdf("figures/strata_by_apoe/brain_mets_strata_class_legend.pdf", height = 5,width = 5)
grid.draw(p_class_legend)
dev.off()

blank <- ggplot()+theme_void()
pdf("figures/strata_by_apoe/brain_mets_tangle_strata_sig.pdf",
    height = 4.5,width = 10,onefile = F) 
egg::ggarrange(
  p_fp,blank,p_class,
  ncol=1, heights =  c(10,-0.3,0.4)) 
dev.off()


## Scatter plot all-----
brain_targeted_info
dat_brain$apoe4_b2 <- as.factor(dat_brain$apoe4_b2)
dat_brain_nona <- dat_brain %>% 
  filter(!is.na(apoe4_b2))

head(brain_targeted_info)
(scatter_list <- unique(re_tang_comb$com_id))

for (i in scatter_list) {
  met_name <- brain_targeted_info$Compounds[which(brain_targeted_info$com_id==i)]
  p_scatter <- 
    ggplot(dat_brain_nona,aes(x=dat_brain_nona[[i]],y=tangsqrt,color=apoe4_b2))+
    geom_point(alpha=0.5)+
    geom_smooth(method = "lm")+
    scale_color_manual(values = c("#338ec4","#c96353"),
                       labels = c("APOE4 non-carriers","APOE4 carriers"),
                       name="")+
    theme_classic()+
    theme(axis.text = element_text(size = 15),
          axis.title = element_blank(),
          legend.text = element_text(size = 18),
          legend.position = "none",
          legend.justification = c(0, 1))
  pdf(paste0("figures/scatter_by_apoe/scatter_tangle_",met_name,".pdf"),
      height = 3.5,width = 4)
  print(p_scatter)
  dev.off()
}




