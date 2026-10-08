rm(list = ls())
library(tidyverse)
library(data.table)
library(openxlsx)
library(sjmisc)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(circlize)
library(ComplexUpset)
library(RColorBrewer)
library(ggsci)
library(scales)
library(ggrepel)



setwd("workpath")
select <- dplyr::select

load("data/dat_combine_all.RData")
re_brain_all <- read.xlsx("results/brain_mets_outcome/brain_mets_outcome.xlsx")

# Volcano plot outcomes -------
re_brain_all <- re_brain_all %>% 
  mutate(log10pvalue= -log10(P),
         labels=case_when(
           fdr<0.05 ~ Compounds,
           TRUE ~ NA_character_),
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
table(re_brain_all$significance)
re_brain_all$significance <- 
  factor(re_brain_all$significance,
         levels = c("FDR<0.05&Beta>0",
                    "P<0.05&FDR>=0.05&Beta>0",
                    "FDR<0.05&Beta<0",
                    "P<0.05&FDR>=0.05&Beta<0",
                    "Non-significant"))
table(re_brain_all$outcome)

(outcome_var <- c("amylsqrt","tangsqrt"))
(outcome_label <- c("Amyloid-β burden","Tau tangle density"))

scales::show_col(c("#9C0824","#D94602","#26456E","#3482B6","grey"))
volcano_color <- 
  c("Positive"="#9C0824",
    "Negative"="#26456E",
    "Non-significant"="grey")

i=1
for (i in seq_along(outcome_var)){
  dat_vol <- re_brain_all %>% 
    filter(outcome== outcome_var[i] )
  (max_abs_value=max(abs(dat_vol$Beta)))
  (fdr_threshold <- max(dat_vol$P[dat_vol$fdr <= 0.05]))
  pdf(paste0("figures/volcano_plot/volcano_brain_",outcome_var[i],".pdf"),
      width = 7,height = 5)
  p1 <- ggplot(dat_vol,aes(x=Beta,y=log10pvalue))+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(fdr_threshold),
               linetype = "dashed", color = "blue") +
    geom_point(aes(color = significance_fdr), alpha = 0.6, size = 3.5)+
    scale_color_manual(values = volcano_color) +
    scale_x_continuous(limits = c(-max_abs_value, max_abs_value))+
    geom_text_repel(aes(label = labels),
                    size = 7,
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.5, "lines"),
                    segment.size = 0.2,
                    segment.color = "black",
                    max.overlaps = 8,
                    show.legend = FALSE)+
    labs(x = "Beta coefficients", y = "-Log10(P)") + 
    theme_bw()+
    theme(text=element_text(size = 22),
          axis.title = element_text(size=25),
          panel.grid = element_blank(),
          plot.margin = margin(t = 5.5, r = 10, b = 5.5, l = 5.5, unit = "pt"),
          legend.position = "none",
          legend.title = element_blank())
  print(p1)
  dev.off()
}

p1_legend <- p1+
  theme(legend.position = "right")
p1_legend

p1_lgd <- ggpubr::get_legend(p1_legend)

pdf("figures/volcano_plot/volcano_outcome_legend.pdf", height = 1,width = 2.5)
grid.draw(p1_lgd)
dev.off()

# Results filter ------
outcome1 <- c("tangsqrt","amylsqrt") 
head(re_brain_all);dim(re_brain_all)
length(unique(re_brain_all$com_id)) 
re_brain_sig <- re_brain_all %>% 
  filter(outcome %in% outcome1) %>%
  filter(fdr<0.05)
length(unique(re_brain_sig$com_id)) 
table(re_brain_sig$Class)
write.xlsx(re_brain_sig,
           "results/brain_mets_outcome/ad_pathology_sig_metabolites.xlsx")

re_brain_sig_dedup <- re_brain_sig %>% 
  filter(!duplicated(com_id))
table(re_brain_sig_dedup$Class)

(ad_path_mets_list <- unique(re_brain_sig$com_id))
write_lines(ad_path_mets_list, "lists/ad_path_mets_list.txt")


re_brain_ad_met <- re_brain_all %>% 
  filter(outcome %in% c("amylsqrt", "tangsqrt", "cogng_demog_slope")) %>% 
  filter(com_id %in% unique(re_brain_sig$com_id)) %>% 
  group_by(outcome) %>% 
  mutate(
    fdr = case_when(
      outcome %in% c("amylsqrt", "tangsqrt") ~ fdr,
      outcome == "cogng_demog_slope" ~ p.adjust(P, method = "BH")
    )
  ) %>% 
  ungroup()
head(re_brain_ad_met)

re_brain_cog <- re_brain_ad_met %>% 
  filter(outcome=="cogng_demog_slope")
table(re_brain_cog$fdr<0.05)


# ComplexHeatmap vertical ------
library(circlize)
dat_beta <- re_brain_ad_met %>% 
  select(Compounds,outcome,Beta) %>% 
  pivot_wider(values_from = Beta,names_from = outcome) %>% 
  column_to_rownames(var = "Compounds") 
head(dat_beta);dim(dat_beta) 

summary(dat_beta$amylsqrt)
summary(dat_beta$tangsqrt)
summary(dat_beta$cogng_demog_slope)

dat_p <- re_brain_ad_met %>% 
  select("Compounds","outcome","fdr") %>% 
  pivot_wider(values_from = fdr,names_from = outcome) %>% 
  column_to_rownames(var = "Compounds") 
head(dat_p);dim(dat_p) # 3/63

dat_sig <- re_brain_ad_met %>% 
  mutate(sig.sign=case_when(
    fdr <0.05 ~ "#",
    fdr>=0.05 & P<0.05 ~ "*",
    P>=0.05 ~ ""
  )) %>% 
  select("Compounds","outcome","sig.sign") %>% 
  pivot_wider(values_from = sig.sign,names_from = outcome) %>% 
  column_to_rownames(var = "Compounds") 
head(dat_sig);dim(dat_sig)


cbind(colnames(dat_beta),colnames(dat_p))
cbind(rownames(dat_beta),rownames(dat_p))

table(abs(dat_beta)>0.3)

sig_mets_info <- brain_targeted_info %>% 
  filter(com_id %in% unique(re_brain_sig$com_id))
identical(rownames(dat_beta),sig_mets_info$Compounds) # T

(col_pal=pal_nejm("default",alpha=0.8)(8))
scales::show_col(col_pal)
class_anno = rowAnnotation(
  Class=sig_mets_info$Class,
  col = list(Class = c("Acylcarnitines" = col_pal[1],
                       "Aminoacids" = col_pal[2],
                       "Biogenic Amines" = col_pal[3],
                       "Glycerides" = col_pal[4],
                       "Glycerophospholipids" = col_pal[5],
                       "Sphingolipids" = col_pal[6],
                       # "Cholesterol Esters" = col_pal[7],
                       "Sugars" = col_pal[8])),
  annotation_name_gp = gpar(fontsize=18),
  annotation_name_rot = 45)
draw(class_anno)

red_blue <- rev(brewer.pal(n=11,name="RdBu"))
breaklist <- seq(-0.3,0.3,by=0.001)
col_pals = colorRampPalette(red_blue)(length(breaklist))
col_fun1 <- colorRamp2(breaklist, col_pals)

breaklist2 <- seq(-0.03,0.03,by=0.0001)
col_pals2 <- colorRampPalette(red_blue)(length(breaklist2))
col_fun2 <- colorRamp2(breaklist2, col_pals2)

dat_beta1 <- dat_beta[,-3]
table(abs(dat_beta1)>0.3)
dat_beta1[dat_beta1 > 0.3] <- 0.3
dat_beta1[dat_beta1 < -0.3] <- -0.3

dat_beta2 <- dat_beta[,3,drop=FALSE]
table(dat_beta2 > 0.03)
table(dat_beta2 < -0.03)
dat_beta2[dat_beta2 > 0.03] = 0.03
dat_beta2[dat_beta2 < -0.03] = -0.03


dat_p1 <- dat_p[,-3]
dat_p2 <- dat_p[,3,drop=FALSE]

dat_sig1 <- dat_sig[,-3]
dat_sig2 <- dat_sig[,3,drop=FALSE]

cell_size <- unit(6, "mm")

ht1 = Heatmap(dat_beta1,name = "Beta for AD Pathology",
              col = col_fun1,
              column_names_rot = 45,
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              clustering_method_columns = "ward.D2",
              width = ncol(dat_beta1) * cell_size * 3,
              height = nrow(dat_beta1) * cell_size * 1.2,
              column_labels = c("Amyloid-β","Tau tangles"),
              row_names_gp = gpar(fontsize=18),
              column_names_gp = gpar(fontsize=20),
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (dat_sig1[i,j] == "#") {
                  grid.text("#",x,y,
                            gp=gpar(fontsize = 15,col="black"))
                } else if (dat_sig1[i,j] == "*") {
                  grid.text("*",x,y-unit(0.2/nrow(dat_beta1),"npc"),
                            gp=gpar(fontsize = 22,col="black"))
                } 
                grid.rect(x, y, width, height, 
                          gp = gpar(col = "grey60", fill = NA, lwd = 1))})
ht1

ht2 = Heatmap(dat_beta2,name = "Beta for Cognitive decline",
              col = col_fun2,
              column_names_rot = 45,
              width = ncol(dat_beta2) * cell_size * 3,
              height = nrow(dat_beta2) * cell_size * 1.2,
              right_annotation = class_anno,
              column_labels = c("Cognitive decline"),
              row_names_gp = gpar(fontsize=18),
              column_names_gp = gpar(fontsize=20),
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (dat_sig2[i,j] == "#") {
                  grid.text("#",x,y,
                            gp=gpar(fontsize = 15,col="black"))
                } else if (dat_sig2[i,j] == "*") {
                  grid.text("*",x,y-unit(0.2/nrow(dat_beta2),"npc"),
                            gp=gpar(fontsize = 22,col="black"))
                }
                grid.rect(x, y, width, height, 
                          gp = gpar(col = "grey60", fill = NA, lwd = 1))})
ht2

ht_list = ht1 + ht2
ht_list

pdf("figures/heatmap_cognitive_ad_pathology_vertical.pdf",
    width = 6,height = 20)
draw(ht_list,
     show_heatmap_legend = FALSE)
dev.off()

lgd1 = Legend(col_fun = col_fun1, 
              at = c(-0.3,-0.15,0,0.15,0.3),
              title = "Beta for \nAD pathology",
              title_gp = gpar(fontsize=8,fontface="bold"),
              labels_gp = gpar(col="black",fontsize=7),
              title_position = "topleft")
lgd2 = Legend(col_fun = col_fun2,
              at = c(-0.03,-0.015,0,0.015,0.03),
              title = "Beta for \ncognitive decline",
              title_gp = gpar(fontsize=8,fontface="bold"),
              labels_gp = gpar(col="black",fontsize=7),
              title_position = "topleft")

Class = c("Acylcarnitines" = col_pal[1],
          "Aminoacids" = col_pal[2],
          "Biogenic Amines" = col_pal[3],
          # "Glycerides" = col_pal[4],
          "Glycerophospholipids" = col_pal[5],
          "Sphingolipids" = col_pal[6],
          # "Cholesterol Esters" = col_pal[7],
          "Sugars" = col_pal[8])

class_legend <- Legend(
  title = "Class",
  at = names(Class),
  title_gp = gpar(fontsize=8,fontface="bold"),
  labels_gp = gpar(col="black",fontsize=7),
  legend_gp = gpar(fill = Class))

lgd_list <- packLegend(lgd1, lgd2, class_legend, direction = "vertical")

pdf("figures/heatmap_cognitive_ad_pathology_legend.pdf",
    width = 2,height = 5)
draw(lgd_list)
dev.off()

# ComplexUpset ------
library(ComplexUpset)

amyloid <- re_brain_ad_met %>% 
  filter(outcome=="amylsqrt"&fdr<0.05) %>% pull(com_id)
tau <- re_brain_ad_met %>% 
  filter(outcome=="tangsqrt"&fdr<0.05) %>% pull(com_id)
cog  <- re_brain_ad_met %>% 
  filter(outcome=="cogng_demog_slope" & fdr<0.05) %>% pull(com_id)
all_met <- unique(c(amyloid, tau, cog))
df <- tibble(
  metabolite = all_met,
  `Amyloid-β`= as.integer(all_met %in% amyloid),
  `Tau tangles`  = as.integer(all_met %in% tau),
  `Cogntive decline` = as.integer(all_met %in% cog)
)
head(df)

set_colors <- c(
  `Amyloid-β` = "#CF7B6B",
  `Tau tangles` = "#6C9FC3",
  `Cogntive decline` = "#86AAA3"
)

p_upset <- 
  ComplexUpset::upset(
  df,
  intersect = c("Amyloid-β", "Tau tangles", "Cogntive decline"),
  sort_intersections_by = "degree",
  # sort_intersections = "ascending",
  sort_sets = "ascending",
  set_sizes = upset_set_size(position = 'right'),
  name = "Intersection size",
  height_ratio = 0.4,
  width_ratio = 0.3,
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 5),
      fill = "grey30"
    )
  ),
  queries = list(
    upset_query(set = "Amyloid-β", fill = "#CF7B6B"),
    upset_query(set = "Tau tangles", fill = "#6C9FC3"), 
    upset_query(set = "Cogntive decline", fill = "#86AAA3") 
  ),
  themes = upset_modify_themes(
    list(
      'intersections_matrix' = theme(
        axis.text.y = element_text(size = 15,colour = "black"),
        axis.title.x = element_blank(),
        panel.grid = element_blank()
      ),
      'Intersection size' = theme(
        axis.text.y = element_text(size = 13,color="black"),
        axis.title.y = element_text(size = 15,color = "black"),
        axis.line.x.bottom = element_line(linewidth = 0.6),
        axis.line.y.left = element_line(linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ),
      'sets' = theme(
        axis.text.x = element_text(size = 13,color = "black"),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 15,color = "black"),
        panel.grid = element_blank()
      )
    )
  )
) +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )
print(p_upset)

pdf("figures/upset_amyloid_tau_cog.pdf",height = 4,width = 7)
print(p_upset)
dev.off()

