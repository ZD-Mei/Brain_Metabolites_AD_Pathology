rm(list = ls())
library(tidyverse)
library(data.table)
library(openxlsx)
library(sjmisc)

setwd("workpath")

select <- dplyr::select

load("data/dat_combine_all.RData")

dat_all$ser_visit
summary(dat_all$age_at_visit)
summary(dat_all$age_death)

# serum mets and outcomes ------
outcome_con <- 
  c("amylsqrt","tangsqrt","cogng_demog_slope") 
adj_vars <- 
  c("cogng_demog_slope")

(i=serum_targeted_list[2])
(j=outcome_con[2])
re_serum_con <- c()
for (i in serum_targeted_list) {
  for (j in outcome_con) {
    if (j %in% adj_vars) {
      mod.form <- as.formula(paste0(j,"~",i,"+age_at_visit+age_death+bmi+phys5itemsum")) 
    } else {
      mod.form <- as.formula(paste0(j,"~",i,"+age_at_visit+age_death+msex+educ+bmi+phys5itemsum")) 
    }
    mod.log <- lm(mod.form,data = dat_all)
    re_serum_con <- rbind(re_serum_con,c(i,j,coef(summary(mod.log))[2,c(1,2,4)]))
  }
}
re_serum_con <- as.data.frame(re_serum_con)
head(re_serum_con);dim(re_serum_con) 
colnames(re_serum_con) <- c("com_id","outcome","Beta","SE","P")
re_serum_con[,3:5] <- map_df(re_serum_con[,3:5],as.numeric)

re_serum_all <- re_serum_con %>% 
  inner_join(serum_targeted_info,by="com_id") %>% 
  group_by(outcome) %>% 
  mutate(fdr=p.adjust(P, n = length(P),method = "BH")) %>% 
  ungroup(.)
head(re_serum_all) 
write.xlsx(re_serum_all,"results/serum_mets_outcome/serum_mets_outcome.xlsx")

# output STable --------
re_serum_stable <- re_serum_all %>% 
  select(-c(com_id,HMDB_main)) %>% 
  select(Compounds,Class,everything())
re_serum_stable_wide <- re_serum_stable %>% 
  pivot_wider(values_from = c(Beta,SE,P,fdr),
              names_from = outcome,
              names_vary = "slowest")
head(re_serum_stable_wide)
write.xlsx(re_serum_stable_wide,
           "results/serum_mets_outcome/stable_serum_mets_outcome.xlsx")

# Volcano plot -------
head(re_serum_all);dim(re_serum_all)
re_serum_all <- re_serum_all %>% 
  mutate(log10pvalue= -log10(P),
         labels=case_when(
           P<0.05 ~ Compounds,
           TRUE ~ NA_character_),
         label_overlap_brain = case_when(
           (Compounds %in% brain_targeted_info$Compounds)&P<0.05 ~ Compounds,
           TRUE ~ NA_character_),
         overlap_brain = case_when(
           Compounds %in% brain_targeted_info$Compounds ~ "Yes",
           TRUE ~ "No"),
         class_overlap_brain = case_when(
           Compounds %in% brain_targeted_info$Compounds ~ Class,
           TRUE ~ "Serum only"),
         significance=case_when(
           fdr<0.05 & Beta >0 ~ "FDR<0.05&Beta>0",
           fdr<0.05 & Beta <0 ~ "FDR<0.05&Beta<0",
           P<0.05 & fdr>=0.05 & Beta >0 ~ "P<0.05&FDR>=0.05&Beta>0",
           P<0.05 & fdr>=0.05 & Beta <0 ~ "P<0.05&FDR>=0.05&Beta<0",
           TRUE ~ "Non-significant"),
         significance_fdr=case_when(
           fdr<0.05 & Beta >0 ~ "Positive",
           fdr<0.05 & Beta <0 ~ "Negative",
           TRUE ~ "Non-significant"),
  )
table(re_serum_all$significance)
re_serum_all$significance <- 
  factor(re_serum_all$significance,
         levels = c("FDR<0.05&Beta>0",
                    "P<0.05&FDR>=0.05&Beta>0",
                    "FDR<0.05&Beta<0",
                    "P<0.05&FDR>=0.05&Beta<0",
                    "Non-significant"))
table(re_serum_all$outcome)

(outcome_var <- unique(re_serum_all$outcome))
(outcome_label <- c("Cognitive decline","Amyloid-beta","Tau tangles"))

volcano_color <- c(
  "FDR<0.05&Beta<0" = "#26456E",
  "P<0.05&FDR>=0.05&Beta<0" = "#3482B6",
  "FDR<0.05&Beta>0" = "#9C0824",
  "P<0.05&FDR>=0.05&Beta>0" = "#D94602",
  "Non-significant" = "grey"
  )


i=1
for (i in seq_along(outcome_var)){
  dat_vol <- re_serum_all %>% 
    filter(outcome== outcome_var[i] )
  (max_abs_value=max(abs(dat_vol$Beta)))
  (fdr_threshold <- max(dat_vol$P[dat_vol$fdr <= 0.05]))
  pdf(paste0("figures/serum_met_ad/volcano_serum_",outcome_var[i],".pdf"),
      width = 4,height = 4)
  p1 <- ggplot(dat_vol,aes(x=Beta,y=log10pvalue))+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(fdr_threshold),
               linetype = "dashed", color = "blue") +
    geom_point(aes(color = significance), alpha = 0.6, size = 1.5)+
    # scale_color_manual(values = c("#9C0824","#D94602","#26456E","#3482B6","grey")) +
    scale_color_manual(values = volcano_color) +
    scale_x_continuous(limits = c(-max_abs_value, max_abs_value))+
    geom_text_repel(aes(label = labels),
                     size = 5,
                     box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.5, "lines"),
                     segment.size = 0.2,
                     segment.color = "black",
                     max.overlaps = 8,
                     show.legend = FALSE)+
    labs(x = "Beta coefficients", y = "-Log10(P)") + #, title = outcome_label[i]
    theme_bw()+
    theme(text=element_text(size = 18),
          # title = element_text(size = 20),
          axis.title = element_text(size=18),
          panel.grid = element_blank(),
          plot.margin = margin(t = 5.5, r = 10, b = 5.5, l = 5.5, unit = "pt"),
          legend.position = "none",
          legend.title = element_blank())
  print(p1)
  dev.off()
}


library(ComplexHeatmap)
lgd <- Legend(
  type = "grid", #"points"
  # pch  = 20, 
  at = c("Negative (FDR<0.05)", "Negative (P<0.05&FDR>0.05)", 
         "Non-significant", "Positive (P<0.05&FDR>0.05)","Positive (FDR<0.05)"),
  legend_gp = gpar(fill = c("#26456E","#3482B6","grey","#D94602","#9C0824"),
                   col = NA
                   ),
  labels_gp = gpar(fontsize=10),
  grid_width = unit(4, "mm"), 
  grid_height = unit(4, "mm"), 
  background = NULL
)

pdf("figures/serum_met_ad/volcano_outcome_legend.pdf",
    height = 1.5,width = 3)
draw(lgd,
     x = unit(0.5, "npc"),
     y = unit(0.5, "npc"),
     just = "center")
dev.off()


