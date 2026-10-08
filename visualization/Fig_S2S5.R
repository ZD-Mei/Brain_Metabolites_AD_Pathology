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
library(RColorBrewer)
library(ggsci)
library(scales)
library(ggrepel)
library(paletteer)


setwd("workpath")
select <- dplyr::select

load("data/dat_combine_all.RData")


# Fig S2 -----
(colnames(dat_brain))

df_brain <- dat_brain[,c(1,78:254)] %>% 
  column_to_rownames(var = "projid")
cor_mat <- cor(df_brain, use = "pairwise.complete.obs", method = "pearson") 

cor_df <- as.data.frame(cor_mat)
cor_df[cor_df==1] <- 0
min(cor_df)
max(cor_df)

(col_pal=pal_nejm("default",alpha=0.8)(8))
scales::show_col(col_pal)
table(brain_targeted_info$Class)
identical(brain_targeted_info$com_id,colnames(cor_mat))
colnames(cor_mat) <- brain_targeted_info$Compounds
rownames(cor_mat) <- brain_targeted_info$Compounds

class_anno = columnAnnotation(
  Class=brain_targeted_info$Class,
  col = list(Class = c("Acylcarnitines" = col_pal[1],
                       "Aminoacids" = col_pal[2],
                       "Biogenic Amines" = col_pal[3],
                       "Glycerides" = col_pal[4],
                       "Glycerophospholipids" = col_pal[5],
                       "Sphingolipids" = col_pal[6],
                       # "Cholesterol Esters" = col_pal[7],
                       "Sugars" = col_pal[8])))

## heatmap with names ------
ht_brain_v2 <- Heatmap(cor_mat, 
                       name = "Correlation", 
                       bottom_annotation = class_anno,
                       cluster_rows = TRUE, 
                       cluster_columns = TRUE)

pdf("figures/cor_heatmap/heatmap_brain_mets_v2.pdf",
    width = 25,height = 25)
draw(ht_brain_v2,
     show_annotation_legend = FALSE,
     show_heatmap_legend = FALSE
)
dev.off()

## output correlations dataframe -----
row_idx <- row_order(ht_brain)
col_idx <- column_order(ht_brain)

row_names_ordered <- rownames(cor_mat)[row_idx]
col_names_ordered <- colnames(cor_mat)[col_idx]

cor_mat <- cor_mat[row_names_ordered,col_names_ordered]
cor_df <- as.data.frame(cor_mat) %>% 
  rownames_to_column(var = "Metabolites")
head(cor_df[,1:5])
write.xlsx(cor_df,"results/basic_info/correlations_brain_mets.xlsx")

# Fig. S5 -----
(colnames(dat_serum))
df_serum <- dat_serum[,c(1,255:524)] %>% 
  column_to_rownames(var = "projid")
cor_mat <- cor(df_serum, use = "pairwise.complete.obs", method = "pearson") 
# cor_mat[upper.tri(cor_mat)] <- NA

(col_pal=pal_nejm("default",alpha=0.8)(8))
scales::show_col(col_pal)
table(serum_targeted_info$Class)
identical(serum_targeted_info$com_id,colnames(cor_mat))
colnames(cor_mat) <- serum_targeted_info$Compounds
rownames(cor_mat) <- serum_targeted_info$Compounds

class_anno = columnAnnotation(
  Class=serum_targeted_info$Class,
  col = list(Class = c("Acylcarnitines" = col_pal[1],
                       "Aminoacids" = col_pal[2],
                       "Biogenic Amines" = col_pal[3],
                       "Glycerides" = col_pal[4],
                       "Glycerophospholipids" = col_pal[5],
                       "Sphingolipids" = col_pal[6],
                       "Cholesterol Esters" = col_pal[7],
                       "Sugars" = col_pal[8])))

## heatmap with names ------
ht_serum_v2 <- Heatmap(cor_mat, 
                       name = "Correlation", 
                       # show_column_names = FALSE,
                       # show_row_names = FALSE,
                       bottom_annotation = class_anno,
                       cluster_rows = TRUE, 
                       cluster_columns = TRUE)

pdf("figures/cor_heatmap/heatmap_serum_mets_v2.pdf",
    width = 30,height = 30)
draw(ht_serum_v2,
     show_annotation_legend = FALSE,
     show_heatmap_legend = FALSE
)
dev.off()

## output correlations dataframe -----
row_idx <- row_order(ht_serum)
col_idx <- column_order(ht_serum)

# Step 4
row_names_ordered <- rownames(cor_mat)[row_idx]
col_names_ordered <- colnames(cor_mat)[col_idx]

cor_mat <- cor_mat[row_names_ordered,col_names_ordered]
cor_df <- as.data.frame(cor_mat) %>% 
  rownames_to_column(var = "Metabolites")
head(cor_df[,1:5])
write.xlsx(cor_df,"results/basic_info/correlations_serum_mets.xlsx")







