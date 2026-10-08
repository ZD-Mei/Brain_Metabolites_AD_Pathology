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

# association results -----
re_assoc_all <- read.xlsx("results/association/assoc_brain_mets_AD.xlsx")
table(re_assoc_all$outcome)
re_brain_assoc <- re_assoc_all %>% 
  filter(outcome %in% c("amylsqrt","tangsqrt","cogng_demog_slope"))

re_amy_sig <- re_assoc_all %>% 
  filter(outcome == "amylsqrt" & fdr <0.05)
dim(re_amy_sig)

re_tang_sig <- re_assoc_all %>% 
  filter(outcome == "tangsqrt" & fdr <0.05)
dim(re_tang_sig)

# apoe4 interaction -------
outcome_con <-  c("amylsqrt","tangsqrt","cogng_demog_slope")
adj_vars <- 
  c("cogng_demog_slope")

(i=brain_targeted_list[2])
(j=outcome_con[2])
re_apoe4_int_list <- list()
k=1
for (i in brain_targeted_list) {
  for (j in outcome_con) {
    if (i %in% adj_vars){
      mod.form <- as.formula(paste0(j,"~",i,"*apoe4_b2+age_death+pmi")) 
    } else {
      mod.form <- as.formula(paste0(j,"~",i,"*apoe4_b2+age_death+msex+educ+pmi")) 
    }
    mod <- lm(mod.form,data = dat_brain)
    summary_mod <- summary(mod)
    var_exp <- c(i,"apoe4_b2",paste0(i,":","apoe4_b2"))
    re_apoe4_int_list[[k]] <- tibble(
      com_id = i,
      outcome = j,
      vars = c("met","apoe4","int_term"),
      Beta = summary_mod$coefficients[var_exp,"Estimate"],
      SE = summary_mod$coefficients[var_exp,"Std. Error"],
      P = summary_mod$coefficients[var_exp,"Pr(>|t|)"],
    )
    k = k+1
  }
}
re_apoe4_int <- bind_rows(re_apoe4_int_list)
head(re_apoe4_int);dim(re_apoe4_int) #1593

re_apoe4_int <- re_apoe4_int %>% 
  filter(vars=="int_term") %>% 
  group_by(outcome) %>% 
  mutate(fdr=p.adjust(P, n = length(P),method = "BH")) %>% 
  ungroup(.) %>% 
  inner_join(brain_targeted_info,by="com_id")
head(re_apoe4_int) 
write.xlsx(re_apoe4_int,
           "results/interaction_apoe/brain_mets_apoe4_int_AD.xlsx")

# Interaction results -----
re_apoe4_int <- read.xlsx("results/interaction_apoe/brain_mets_apoe4_int_AD.xlsx")
head(re_apoe4_int);dim(re_apoe4_int)
re_apoe4_int <- re_apoe4_int %>% 
  select(-c(vars,Compounds,Class,HMDB_main))
colnames(re_apoe4_int)[3:6] <- paste0(colnames(re_apoe4_int[3:6]),"_int")

re_apoe4_int_sig <- re_apoe4_int %>% 
  filter(fdr_int <0.05) %>% 
  mutate(pair1=paste0(com_id,"_",outcome))
length(unique(re_apoe4_int_sig$com_id)) #14
re_apoe4_int_sig$pair1

re_int_amy <- re_apoe4_int %>% 
  filter(outcome=="amylsqrt") %>% 
  filter(com_id %in% re_amy_sig$com_id) %>% 
  mutate(fdr_int_amy=p.adjust(P_int,method = "BH")) %>% 
  arrange(fdr_int_amy)
re_int_amy

re_int_tang <- re_apoe4_int %>% 
  filter(outcome=="tangsqrt") %>% 
  filter(com_id %in% re_tang_sig$com_id) %>% 
  mutate(fdr_int_tang=p.adjust(P_int,method = "BH")) %>% 
  arrange(fdr_int_tang)
re_int_tang

# Volcano plot amyloid -----
head(re_int_amy)
re_int_amy_viz <- re_int_amy %>% 
  mutate(lgP_apoe4= -log10(P_int)) %>% 
  left_join(brain_targeted_info,by="com_id") %>% 
  mutate(
    lab_apoe4=case_when(
      P_int<0.05 ~ Compounds,
      TRUE ~ NA_character_),
    sig_apoe4=case_when(
      fdr_int_amy<0.05 & Beta_int >0 ~ "FDR<0.05&Interaction term>0",
      fdr_int_amy<0.05 & Beta_int <0 ~ "FDR<0.05&Interaction term<0",
      P_int<0.05 & fdr_int_amy>=0.05 & Beta_int >0 ~ "P<0.05&FDR>=0.05&Interaction term>0",
      P_int<0.05 & fdr_int_amy>=0.05 & Beta_int <0 ~ "P<0.05&FDR>=0.05&Interaction term<0",
      TRUE ~ "Non-significant"))

re_int_amy_viz$sig_apoe4 <- 
  factor(re_int_amy_viz$sig_apoe4,
         levels = c(
           "FDR<0.05&Interaction term<0",
           "P<0.05&FDR>=0.05&Interaction term<0",
           "P<0.05&FDR>=0.05&Interaction term>0",
           "FDR<0.05&Interaction term>0",
           "Non-significant"))


(max_abs_value=max(abs(re_int_amy_viz$Beta_int)))

col_panel <- c(
  `FDR<0.05&Interaction term<0`="#193b1b",
  `P<0.05&FDR>=0.05&Interaction term<0`="#74b577",
  `P<0.05&FDR>=0.05&Interaction term>0`="#e491b6",
  `FDR<0.05&Interaction term>0`="#602f5a", #"#38123d"
  `Non-significant`="grey50")
show_col(col_panel)

pdf("figures/volcano_plot/volcano_int_apoe4_amyloid_filter.pdf",
    width = 5,height = 3)
p1 <- ggplot(re_int_amy_viz,aes(x=Beta_int,y=lgP_apoe4))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_point(aes(color = sig_apoe4), alpha = 0.8, size = 3)+
  scale_color_manual(values = col_panel) +
  scale_x_continuous(limits = c(-max_abs_value, max_abs_value))+
  geom_text_repel(aes(label = lab_apoe4),
                  size = 6,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.size = 0.2,
                  segment.color = "black",
                  max.overlaps = 12,
                  show.legend = FALSE)+
  theme_bw()+
  theme(text=element_text(size = 18),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 5.5, r = 10, b = 5.5, l = 5.5, unit = "pt"),
        legend.position = "none",
        legend.title = element_blank())
print(p1)
dev.off()

# Volcano plot tangle -----
re_int_tang_viz <- re_int_tang %>% 
  mutate(lgP_apoe4= -log10(P_int)) %>% 
  left_join(brain_targeted_info,by="com_id") %>% 
  mutate(
    lab_apoe4=case_when(
      fdr_int_tang<0.05 ~ Compounds,
      TRUE ~ NA_character_),
    sig_apoe4=case_when(
      fdr_int_tang<0.05 & Beta_int >0 ~ "FDR<0.05&Interaction term>0",
      fdr_int_tang<0.05 & Beta_int <0 ~ "FDR<0.05&Interaction term<0",
      P_int<0.05 & fdr_int_tang>=0.05 & Beta_int >0 ~ "P<0.05&FDR>=0.05&Interaction term>0",
      P_int<0.05 & fdr_int_tang>=0.05 & Beta_int <0 ~ "P<0.05&FDR>=0.05&Interaction term<0",
      TRUE ~ "Non-significant"))

re_int_tang_viz$sig_apoe4 <- 
  factor(re_int_tang_viz$sig_apoe4,
         levels = c(
           "FDR<0.05&Interaction term<0",
           "P<0.05&FDR>=0.05&Interaction term<0",
           "P<0.05&FDR>=0.05&Interaction term>0",
           "FDR<0.05&Interaction term>0",
           "Non-significant"))

(max_abs_value=max(abs(re_int_tang_viz$Beta_int)))
(fdr_threshold <- max(re_int_tang_viz$P_int[re_int_tang_viz$fdr_int_tang <= 0.05]))
table(re_int_tang_viz$sig_apoe4)


pdf("figures/volcano_plot/volcano_int_apoe4_tangle_filter.pdf",
    width = 5,height = 3)
p1 <- ggplot(re_int_tang_viz,aes(x=Beta_int,y=lgP_apoe4))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(fdr_threshold),
             linetype = "dashed", color = "blue") +
  geom_point(aes(color = sig_apoe4), alpha = 0.8, size = 3)+
  scale_color_manual(values = col_panel) +
  scale_x_continuous(limits = c(-max_abs_value, max_abs_value))+
  geom_text_repel(aes(label = lab_apoe4),
                  size = 5,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.size = 0.2,
                  segment.color = "black",
                  max.overlaps = 12,
                  show.legend = FALSE)+
  theme_bw()+
  theme(text=element_text(size = 18),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 5.5, r = 10, b = 5.5, l = 5.5, unit = "pt"),
        legend.position = "none",
        legend.title = element_blank())
print(p1)
dev.off()

p1_lgd <- p1 + 
  guides(color = guide_legend(override.aes = list(size = 4),ncol = 1))+
  theme(legend.position = "right")
print(p1_lgd)

vol_lgd <- ggpubr::get_legend(p1_lgd)
pdf("figures/volcano_plot/volcano_int_legend.pdf",
    width = 4,height = 2)
grid.draw(vol_lgd)
dev.off()

