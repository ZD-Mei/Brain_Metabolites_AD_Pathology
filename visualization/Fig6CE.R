rm(list = ls())
library(tidyverse)
library(data.table)
library(openxlsx)
library(sjmisc)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(ComplexUpset)
library(RColorBrewer)
library(pheatmap)
library(ggsci)
library(scales)
library(ggrepel)


setwd("workpath")

select <- dplyr::select

# load data -----
load("data/dat_combine_all.RData")
head(brain_targeted_info)
head(serum_targeted_info)
colnames(brain_targeted_info) <- 
  c("brain_com_id","brain_compound","brain_Class","brain_HMDB")
colnames(serum_targeted_info) <- 
  c("serum_com_id","serum_compound","serum_Class","serum_HMDB")

brain_serum_info <- brain_targeted_info %>% 
  inner_join(serum_targeted_info,
             by=c("brain_compound"="serum_compound",
                  "brain_Class"="serum_Class")) 
head(brain_serum_info);dim(brain_serum_info) 


# brain met and ad ------
res_brain_ad <- 
  read.xlsx("results/brain_mets_outcome/brain_mets_outcome.xlsx") 
head(res_brain_ad);dim(res_brain_ad)

res_brain_ad_sig <- res_brain_ad %>% 
  filter(outcome != "cogng_demog_slope") %>% 
  filter(fdr<0.05)
length(unique(res_brain_ad_sig$com_id)) # 63
ad_brain_met_list <- unique(res_brain_ad_sig$com_id)

res_brain_ad_filter <- res_brain_ad %>% 
  filter(com_id %in% unique(res_brain_ad_sig$com_id)) %>% 
  mutate(t = Beta/SE,
         sig = case_when(
           fdr<0.05 ~ "#",
           fdr>=0.05 & P<0.05 ~ "*",
           TRUE ~ ""),
         group = case_when(
           fdr<0.05 & Beta >0 ~ "Pos_sig",
           fdr<0.05 & Beta <0 ~ "Neg_sig",
           fdr>=0.05 & P<0.05 & Beta > 0 ~ "Pos_nomsig",
           fdr>=0.05 & P<0.05 & Beta < 0 ~ "Neg_nomsig",
           TRUE ~ "Non_sig"
         ))
summary(res_brain_ad_filter$t)
table(res_brain_ad_filter$group)

res_brain_annot <- res_brain_ad_filter %>% 
  select(Compounds,outcome,group) %>% 
  pivot_wider(values_from = group,names_from = outcome) %>% 
  column_to_rownames(var = "Compounds") %>% 
  setnames(c("Amyloid","Tau tangles","Cognitive decline"))
head(res_brain_annot)

# serum met and ad ------
res_serum_ad <- 
  read.xlsx("results/serum_mets_outcome/serum_mets_outcome.xlsx") 
head(res_serum_ad);dim(res_serum_ad)

res_serum_ad <- res_serum_ad %>% 
  mutate(t = Beta/SE,
         sig = case_when(
           fdr<0.05 ~ "#",
           fdr>=0.05 & P<0.05 ~ "*",
           TRUE ~ ""),
         group = case_when(
           fdr<0.05 & Beta >0 ~ "Pos_sig",
           fdr<0.05 & Beta <0 ~ "Neg_sig",
           fdr>=0.05 & P<0.05 & Beta > 0 ~ "Pos_nomsig",
           fdr>=0.05 & P<0.05 & Beta < 0 ~ "Neg_nomsig",
           TRUE ~ "Non_sig"
         ))
summary(res_serum_ad$t)
table(res_serum_ad$group)

res_serum_annot <- res_serum_ad %>% 
  select(Compounds,outcome,group) %>% 
  pivot_wider(values_from = group,names_from = outcome) %>% 
  column_to_rownames(var = "Compounds") %>% 
  setnames(c("Amyloid","Tau tangles","Cognitive decline"))
head(res_serum_annot)

res_serum_ad_nominal_sig <- res_serum_ad %>% 
  filter(P<0.05)
ad_serum_nomi_sig_list <- unique(res_serum_ad_nominal_sig$com_id)

# VennDiagram ------
head(brain_targeted_info)
head(serum_targeted_info)

scales::show_col(c("#5b92f8","#f36e65","#20af42"))
table(brain_targeted_info$Compounds %in% serum_targeted_info$Compounds)
table(serum_targeted_info$Compounds %in% brain_targeted_info$Compounds)

library(eulerr)
venn_data <- euler(c(
  `Brain` = 21,         # Only A
  `Serum` = 114,         # Only B
  "Brain&Serum" = 156     # Intersection of A and B
))

pdf("figures/brain_serum/brain_serum_metabolites.pdf",height = 4,width = 4)
plot(venn_data,
     fills = list(fill=c("lightblue", "lightcoral"),alpha =0.5),
     edges = list(col = "black", lwd = 2),
     labels = FALSE,
     quantities = FALSE
)
dev.off()

# Association brain & blood---------
colnames(dat_all)
head(brain_targeted_list)
head(serum_targeted_list)

(x=serum_targeted_list[1])
(y=brain_targeted_list[1])
data = dat_all

fit_mod <- function(x, y, data){
  form <- as.formula(paste0(y,"~",x,"+age_death+age_at_visit+pmi+msex"))
  mod <- try(lm(form,data = data),silent = TRUE)
  
  if (inherits(mod, "try-error")) {
    return(tibble(beta = NA_real_, se = NA_real_, t = NA_real_,
                  p = NA_real_, n = NA_integer_))}
  coefs <- summary(mod)$coefficients
  if (!(x %in% rownames(coefs))) {
    return(tibble(beta = NA_real_, se = NA_real_, t = NA_real_,
                  p = NA_real_, n = NA_integer_))
  }
  tibble(
    beta = unname(coefs[x, "Estimate"]),
    se   = unname(coefs[x, "Std. Error"]),
    t   = unname(coefs[x, "t value"]),
    p    = unname(coefs[x, "Pr(>|t|)"]),
    n    = nobs(mod)
  )
}

comb_index <- tidyr::crossing(
  serum_com_id = serum_targeted_list,
  brain_com_id = brain_targeted_list
)

res_asso_brain_serum <- 
  purrr::pmap_dfr(comb_index, ~ fit_mod(..1, ..2, data = dat_all)) %>% 
  dplyr::bind_cols(comb_index,.) 
head(res_asso_brain_serum)

res_asso_brain_serum2 <- res_asso_brain_serum %>% 
  mutate(fdr=p.adjust(p,method = "BH"),
         sig=case_when(
           fdr<0.05 ~ "#",
           fdr>=0.05 & p<0.05 ~ "*",
           TRUE ~ "")) %>% 
  left_join(brain_targeted_info,by = c("brain_com_id")) %>% 
  left_join(serum_targeted_info,by = c("serum_com_id")) %>%
  arrange(p)
head(res_asso_brain_serum2)
write.xlsx(res_asso_brain_serum2,
           "results/brain_serum/association_brain_serum_mets.xlsx")

res_asso_brain_serum_wide <- res_asso_brain_serum2 %>% 
  filter(brain_com_id %in% ad_brain_met_list) %>% 
  arrange(brain_com_id,serum_com_id) %>% 
  select(-c(serum_com_id,brain_com_id,t,sig,n,brain_Class,brain_HMDB,serum_HMDB)) %>% 
  pivot_wider(values_from = c(beta,se,p,fdr),names_from = brain_compound,
              names_vary = "slowest")

wb <- createWorkbook()
addWorksheet(wb, "sheet1")
writeData(wb, sheet = 1, res_asso_brain_serum_wide, startCol = 1, 
          startRow = 3, rowNames = FALSE, colNames = FALSE)
wb$worksheets[[1]]$mergeCells <- NULL

names <- colnames(res_asso_brain_serum_wide)[3:ncol(res_asso_brain_serum_wide)]
metabolites <- gsub("^(beta|se|p|fdr)_", "", names)
unique_mets <- unique(metabolites)

# col_start <- 1
(met <- unique_mets[1])
for (met in unique_mets) {
  cols <- which(metabolites == met)
  excel_cols <- cols + 2
  mergeCells(wb, sheet=1, rows=1, cols = min(excel_cols):max(excel_cols))     
  writeData(wb, 1, met, startCol=min(excel_cols), startRow=1)
  
  # beta/se/p/fdr
  stat_names <- gsub("_.+$", "", names[cols])
  writeData(wb, 1, t(as.matrix(stat_names)), startCol=min(excel_cols), 
            startRow=2, colNames = FALSE, rowNames = FALSE)
}
saveWorkbook(wb, "results/brain_serum/stable_brain_serum_asso_wide.xlsx",
             overwrite = TRUE)

## heatmap  -----
res_asso_brain_serum_sig <- res_asso_brain_serum2 %>%
  filter(brain_com_id %in% ad_brain_met_list) %>%
  filter(fdr<0.05)
length(unique(res_asso_brain_serum_sig$brain_com_id)) #27
length(unique(res_asso_brain_serum_sig$serum_com_id)) #47
summary(res_asso_brain_serum_sig$beta)

table(unique(res_asso_brain_serum_sig$serum_com_id) %in% ad_serum_nomi_sig_list)

res_asso_wide_sig <- res_asso_brain_serum2 %>% 
  filter(brain_com_id %in% unique(res_asso_brain_serum_sig$brain_com_id)) %>% 
  filter(serum_com_id %in% unique(res_asso_brain_serum_sig$serum_com_id)) %>% 
  select(brain_compound,serum_compound,beta,sig) %>% 
  pivot_wider(values_from = c(beta,sig),names_from = serum_compound) %>% 
  column_to_rownames(var = "brain_compound")
head(res_asso_wide_sig[1:10])
# rows are brain; cols are serum

res_r_sig <- res_asso_wide_sig %>% 
  select(grep("^beta_",colnames(res_asso_wide_sig))) %>% as.matrix(.)
colnames(res_r_sig) <- sub("^beta_","",colnames(res_r_sig))
head(res_r_sig[1:10]);dim(res_r_sig)
colnames(res_r_sig)

res_star_sig <- res_asso_wide_sig %>% 
  select(grep("^sig_",colnames(res_asso_wide_sig))) %>% as.matrix(.)
colnames(res_star_sig) <- sub("^sig_","",colnames(res_star_sig))

## plot  ------
head(rownames(res_r_sig))
head(colnames(res_r_sig))
head(res_brain_annot)
brain_annot <- res_brain_annot[rownames(res_r_sig),]
serum_annot <- res_serum_annot[colnames(res_r_sig),]

legend_ord <- c("Neg_sig","Neg_nomsig","Non_sig","Pos_nomsig","Pos_sig")

row_annot <- rowAnnotation(
  Amyloid = brain_annot$Amyloid,
  `Tau tangles` = brain_annot$`Tau tangles`,
  `Cognitive decline` = brain_annot$`Cognitive decline`,
  col = list(
    Amyloid = c("Neg_sig"="#2E5A87","Neg_nomsig"="#6B9AC2",
                "Pos_nomsig"="#F26A5D","Pos_sig"="#B7203C","Non_sig"="grey"),
    `Tau tangles` = c("Neg_sig"="#2E5A87","Neg_nomsig"="#6B9AC2",
                      "Pos_nomsig"="#F26A5D","Pos_sig"="#B7203C","Non_sig"="grey"),
    `Cognitive decline` = c("Neg_sig"="#2E5A87","Neg_nomsig"="#6B9AC2",
                          "Pos_nomsig"="#F26A5D","Pos_sig"="#B7203C","Non_sig"="grey")
  ),
  annotation_legend_param = list(
    Amyloid = list(at = legend_ord),
    `Tau tangles` = list(at = legend_ord),
    `Cognitive decline` = list(at = legend_ord)),
  annotation_name_gp = gpar(
    fontsize = 16)
  )

col_annot <- HeatmapAnnotation(
  Amyloid = serum_annot$Amyloid,
  `Tau tangles` = serum_annot$`Tau tangles`,
  `Cognitive decline` = serum_annot$`Cognitive decline`,
  annotation_name_side = "left",
  col = list(
    Amyloid = c("Neg_sig"="#2E5A87","Neg_nomsig"="#6B9AC2",
                "Pos_nomsig"="#F26A5D","Pos_sig"="#B7203C","Non_sig"="grey"),
    `Tau tangles` = c("Neg_sig"="#2E5A87","Neg_nomsig"="#6B9AC2",
                      "Pos_nomsig"="#F26A5D","Pos_sig"="#B7203C","Non_sig"="grey"),
    `Cognitive decline` = c("Neg_sig"="#2E5A87","Neg_nomsig"="#6B9AC2",
                          "Pos_nomsig"="#F26A5D","Pos_sig"="#B7203C","Non_sig"="grey")),
  annotation_legend_param = list(
    Amyloid = list(at = legend_ord),
    `Tau tangles` = list(at = legend_ord),
    `Cognitive decline` = list(at = legend_ord)),
  annotation_name_gp = gpar(
    fontsize = 16)
  )


col_fun1 = colorRamp2(seq(-0.3, 0.3, length.out = 7), rev(brewer.pal(7, "RdYlBu")))

ht_v2 <- Heatmap(res_r_sig,
                 name = "Beta",
                 col = col_fun1,
                 row_title = "Brain metabolites",
                 column_title = "Serum metabolites",
                 row_names_gp = gpar(fontsize = 16),
                 column_names_gp = gpar(fontsize = 16),
                 row_title_gp = gpar(fontsize = 18,fontface="bold"),
                 column_title_gp = gpar(fontsize = 18,fontface="bold"),
                 width = ncol(res_r_sig)*unit(5,"mm"),
                 height = nrow(res_r_sig)*unit(5,"mm"),
                 right_annotation = row_annot,
                 bottom_annotation = col_annot,
                 heatmap_legend_param = list(
                   at     = c(-0.3, -0.15, 0, 0.15, 0.3),
                   labels = c("-0.3", "-0.15", "0", "0.15", "0.3")
                 ),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(res_star_sig[i,j],x,
                             ifelse(res_star_sig[i,j]=="*",
                                    y - unit(0.15/nrow(res_star_sig),"npc"),y),
                             gp = gpar(fontsize = 
                                         ifelse(res_star_sig[i,j]=="*",16,10)))},
                 rect_gp = gpar(col="grey70",lwd=0.6)
                 )
print(ht_v2)

pdf("figures/brain_serum/ht_brain_serum_assoc_filtered_annotation.pdf",
    height = 9,width = 15)
draw(ht_v2,show_heatmap_legend = FALSE,show_annotation_legend=FALSE)
dev.off()

# legend 
lgd <- Legend(
  at = c("Negative (FDR<0.05)", "Negative (P<0.05&FDR>0.05)", 
         "Non-significant", "Positive (P<0.05&FDR>0.05)","Positive (FDR<0.05)"),
  legend_gp = gpar(fill = c("#2E5A87","#6B9AC2","grey","#F26A5D","#B7203C"))
)

pdf("figures/brain_serum/legend_annot_brain_serum_assoc.pdf",
    height = 1.5,width = 3)
draw(lgd, 
     x = unit(0.5, "npc"),   
     y = unit(0.5, "npc"),   
     just = "center")
dev.off()

lgd <- Legend(
  title = "Beta",
  col_fun = col_fun1,
  at = seq(-0.3, 0.3, length.out = 5),
  labels = seq(-0.3, 0.3, length.out = 5)
)

pdf("figures/brain_serum/legend_heatmap_brain_serum_assoc.pdf",
    height = 1.5,width = 1)
draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
dev.off()

## 1-on-1 -------
head(brain_serum_info)

i=1
re_assoc <- c()
for (i in seq(nrow(brain_serum_info))) {
  brain_met <- brain_serum_info$brain_com_id[i]
  serum_met <- brain_serum_info$serum_com_id[i]
  (mod_form <- 
      reformulate(termlabels = c(serum_met,"age_death","age_at_visit","pmi","msex"),
                  response = brain_met))
  mod <- lm(mod_form,data = dat_all)
  re_assoc <- rbind(re_assoc,
                    c(brain_serum_info$brain_compound[i],
                      coef(summary(mod))[serum_met,]))
}
head(re_assoc);dim(re_assoc) #156
re_assoc <- as.data.frame(re_assoc)
colnames(re_assoc) <- c("brain_compound","beta","se","t","p")
re_assoc[,2:5] <- map_df(re_assoc[,2:5],as.numeric)
head(re_assoc)
summary(re_assoc$beta)
re_assoc$fdr <- p.adjust(re_assoc$p,method = "BH")
re_assoc <- re_assoc %>% 
  left_join(brain_serum_info,by=c("brain_compound")) %>% 
  mutate(significance=case_when(
    fdr<0.05 ~ "#",
    fdr>=0.05 & p<0.05 ~ "*",
    TRUE ~ NA_character_),
    t2=case_when(
      t >0 ~ t + 0.2,
      t <0 ~ t - 0.2),
    beta2=case_when(
      beta >0 ~ beta + 0.015,
      beta <0 ~ beta - 0.02),
    Class2 = case_when(
      p<0.05 ~ brain_Class,
      TRUE ~ NA_character_
    )) %>% 
  arrange(beta2)

## output STable 
head(re_assoc)
re_assoc_stable <- re_assoc %>% 
  select(brain_compound,brain_Class,beta,se,t,p,fdr) %>% 
  arrange(brain_compound,brain_Class)
head(re_assoc_stable)
write.xlsx(re_assoc_stable,
           "results/brain_serum/stable_brain_serum_association_1on1.xlsx")


### barplot --------
re_assoc$brain_compound <- 
  factor(re_assoc$brain_compound,levels = re_assoc$brain_compound)

(col_pal=pal_nejm("default",alpha=0.8)(8))
scales::show_col(col_pal)
class_col = c("Acylcarnitines" = col_pal[1],
              "Aminoacids" = col_pal[2],
              "Biogenic Amines" = col_pal[3],
              "Glycerides" = col_pal[4],
              "Glycerophospholipids" = col_pal[5],
              "Sphingolipids" = col_pal[6],
              # "Cholesterol Esters" = col_pal[7],
              "Sugars" = col_pal[8])

p_assoc <- ggplot(re_assoc,aes(x=brain_compound,y=beta,fill=Class2))+
  geom_col()+
  # geom_text(aes(y=beta2,label=significance),size=3)+
  scale_fill_manual(values = class_col,na.value = "grey60",name="Class")+
  labs(x="Compounds",
       y="Association of the same metabolites\n measured in brain and serum")+
  theme_classic()+
  theme(text = element_text(size = 15),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank()
  )
print(p_assoc)
ggsave(p_assoc,filename = "figures/brain_serum/association_brain_serum_1on1.pdf",
       width = 8,height = 5)

p_assoc_lgd <- p_assoc + 
  theme(legend.position = "right")+
  guides(fill = guide_legend(nrow = 4))
p_assoc_legend <- ggpubr::get_legend(p_assoc_lgd)

pdf("figures/brain_serum/association_brain_serum_all_legend.pdf",
    width = 7,height = 2)
grid.draw(p_assoc_legend)
dev.off()


## significant barplots -----
re_assoc_sig <- re_assoc %>% 
  filter(p<0.05)
head(re_assoc_sig);dim(re_assoc_sig)

p_assoc_sig <- ggplot(re_assoc_sig,aes(x=brain_compound,y=beta,fill=brain_Class))+
  geom_col()+
  geom_text(aes(y=beta2,label=significance),size=6)+
  scale_fill_manual(values = class_col,na.value = "grey60",name="Class")+
  labs(y="Peason correlation")+
  theme_classic()+
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")
print(p_assoc_sig)
ggsave(p_assoc_sig,
       filename = "figures/brain_serum/association_brain_serum_sig.pdf",
       width = 7,height = 5)

p_assoc_sig_pos <- re_assoc_sig %>% 
  filter(beta > 0) %>% 
  ggplot(aes(x=brain_compound,y=beta,fill=brain_Class))+
  geom_col()+
  geom_text(aes(y=beta2,label=significance),size=8)+
  scale_fill_manual(values = class_col,na.value = "grey60",name="Class")+
  labs(y="Beta coeficient")+
  theme_classic()+
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,
                                   colour = "black",size = 20),
        axis.title.x = element_blank(),
        legend.position = "none")
print(p_assoc_sig_pos)
ggsave(p_assoc_sig_pos,
       filename = "figures/brain_serum/association_brain_serum_sig_positive.pdf",
       width = 7,height = 3.5)

p_assoc_sig_neg <- re_assoc_sig %>% 
  filter(beta < 0) %>% 
  ggplot(aes(x=brain_compound, y=beta, fill= brain_Class))+
  geom_col()+
  geom_text(aes(y=beta2,label=significance),size=8)+
  scale_fill_manual(values = class_col,na.value = "grey60",name="Class")+
  labs(y="Beta coeficient")+
  theme_classic()+
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,
                                   colour = "black",size = 20),
        axis.title.x = element_blank(),
        legend.position = "none")
print(p_assoc_sig_neg)
ggsave(p_assoc_sig_neg,
       filename = "figures/brain_serum/association_brain_serum_sig_negative.pdf",
       width = 3,height = 3.5)




