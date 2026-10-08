rm(list = ls())

library(tidyverse)
library(data.table)
library(openxlsx)
library(coloc)
library(ggsci)
library(paletteer)
library(grid)

select <- dplyr::select

setwd("workpath")

load("data/dat_combine_all.RData")
head(brain_targeted_info);dim(brain_targeted_info)
table(brain_targeted_info$Class)

re_coloc_sig_filter <- 
  read.xlsx("coloc/res_combine/res_eqtl_coloc_sig_filter.xlsx")
head(re_coloc_sig_filter)
table(re_coloc_sig_filter$gene)
re_coloc_sig_filter <- re_coloc_sig_filter %>% 
  filter(!gene %in% c("AC009093.2","AC093298.2","AC119403.1"))

re_coloc_viz <- re_coloc_sig_filter %>% 
  select(Compounds,celltype,gene,PP.H4.abf)

plot_df <- re_coloc_viz %>%
  group_by(Compounds, celltype) %>%
  summarise(
    gene_label = paste(unique(gene), collapse = "\n"),
    n_gene = n_distinct(gene),
    max_pph4 = max(PP.H4.abf,na.rm = TRUE),
    .groups = "drop"
  )
plot_df


# make a full grid so empty cells are shown
plot_df_full <- plot_df %>%
  complete(Compounds, celltype, fill = list(gene_label = "", n_gene = 0)) %>% 
  left_join(brain_targeted_info,by="Compounds")

table(plot_df_full$Class)
plot_df_full$Class <- 
  factor(plot_df_full$Class,
         levels = c("Aminoacids","Acylcarnitines",
                    "Glycerophospholipids","Sphingolipids"))
plot_df_full <- plot_df_full %>% 
  arrange(Compounds) %>% 
  arrange(Class)

# optional: set order of columns
celltype_order <- c("Exc","Inh","Ast","Mic","Oli","End","OPC")
celltype_label <- c("Exc.","Inh.","Ast.","Mic.","Oli.","End.","OPC.")
plot_df_full$celltype <- 
  factor(plot_df_full$celltype, 
         levels = (celltype_order),
         labels = (celltype_label))

# optional: keep compounds in original order or reverse for top-to-bottom display
compound_order <- rev(unique(plot_df_full$Compounds))
plot_df_full$Compounds <- factor(plot_df_full$Compounds, levels = compound_order)

# heatmap
table(plot_df_full$n_gene)
plot_df_full$n_gene <- factor(plot_df_full$n_gene)

head(plot_df_full)
p1 <- ggplot(plot_df_full, aes(x = celltype, y = Compounds, fill = max_pph4)) +
  geom_tile(color = "grey80", linewidth = 0.4) +
  geom_text(aes(label = gene_label), size = 3.5, color = "black",
            fontface = "italic", lineheight = 0.9) +
  scale_fill_viridis_c(option="G",direction=-1,begin=0.6,end = 1,
                       na.value="white",name = "PPH4")+
  labs(x = "Cell types", y = "Metabolites") +
  scale_x_discrete(position = "top")+
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 18,colour = "black"),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.grid = element_blank()
  )
p1

pph4_lgd <- p1+
  theme(legend.position = "right")
pph4_legend <- ggpubr::get_legend(pph4_lgd) 

pdf("figures/coloc/pph4_legend.pdf", height = 2,width = 1)
grid.draw(pph4_legend)
dev.off()


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

p2 <- ggplot(plot_df_full,aes(y=Compounds,x=1,fill = Class))+
  geom_col(width = 1)+
  coord_cartesian(expand = FALSE)+
  scale_fill_manual(values = my_color)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 16,color = "black"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.ticks.x = element_blank())
p2


p_class_lgd <- p2+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
p_class_legend <- ggpubr::get_legend(p_class_lgd) 

pdf("figures/coloc/coloc_met_class_legend.pdf", height = 1,width = 6)
grid.draw(p_class_legend)
dev.off()

blank <- ggplot()+theme_void()
pdf("figures/coloc/coloc_heatmap_horizontal.pdf",
    height = 5,width = 10,onefile = F) 
egg::ggarrange(
  p2,blank,p1,
  nrow=1, widths = c(0.3,-0.3,10)) 
dev.off()


