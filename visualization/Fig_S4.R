rm(list = ls())
library(tidyverse)
library(data.table)
library(openxlsx)
library(sjmisc)
library(grid)
library(gridExtra)
library(missRanger)
library(pheatmap)
library(RColorBrewer)
library(ComplexUpset)
library(RColorBrewer)
library(ggsci)
library(scales)
library(ggrepel)
library(nnet)
library(MASS)
library(CMAverse)
library(paletteer)


setwd("workpath")

select <- dplyr::select

int <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}

# load data -----
load("data/dat_combine_all.RData")

outcome_con <- c("cogng_demog_slope","amylsqrt","tangsqrt")

# _Serum association ------
re_serum <- read.xlsx("results/serum_mets_outcome/serum_mets_outcome.xlsx")
head(re_serum)

# _Brain association ------
re_brain <- read.xlsx("results/brain_mets_outcome/brain_mets_outcome.xlsx")
head(re_brain)

# combine serum & brain ----
head(re_serum)
head(re_brain)

re_serum <- re_serum %>% 
  select(-c(com_id,Class,HMDB_main,fdr)) %>% 
  select(Compounds,outcome,everything())
head(re_serum)
colnames(re_serum)[3:5] <- paste0(colnames(re_serum)[3:5],"_serum")

re_brain <- re_brain %>% 
  select(-c(com_id,fdr)) %>% 
  select(Compounds,Class,HMDB_main,outcome,everything())
head(re_brain)
colnames(re_brain)[5:7] <- paste0(colnames(re_brain)[5:7],"_brain")

## scatter plot -------
re_combine <- re_serum %>% 
  inner_join(re_brain,by=c("Compounds","outcome")) %>% 
  group_by(outcome) %>% 
  mutate(
    fdr_brain=p.adjust(P_brain,method = "BH",n=length(P_brain)),
    fdr_serum=p.adjust(P_serum,method = "BH",n=length(P_serum))) %>% 
  ungroup() %>% 
  mutate(
    Class2=case_when(
      P_serum <0.05 & P_brain <0.05 & Beta_serum * Beta_brain > 0 ~ Class,
      TRUE ~ NA_character_),
    annotation=case_when(
      P_serum<0.05 & P_brain<0.05 & Beta_serum * Beta_brain > 0 ~ Compounds,
      TRUE ~ NA_character_)) %>% 
  select("Compounds","outcome","Beta_serum","SE_serum","P_serum","fdr_serum",
         "Beta_brain","SE_brain","P_brain","fdr_brain",everything())
head(re_combine);dim(re_combine)
length(unique(re_combine$Compounds))
write.xlsx(re_combine,"results/association/assoc_brain_serum_comparison.xlsx")

table(re_combine$outcome)
re_combine_path <- re_combine %>% 
  filter(outcome!="cogng_demog_slope")
head(re_combine_path)
cor.test(re_combine_path$Beta_brain,re_combine_path$Beta_serum)


re_combine_amy <- re_combine %>% 
  filter(outcome=="amylsqrt")
head(re_combine_amy)
cor.test(re_combine_amy$Beta_brain,re_combine_amy$Beta_serum)

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

p_amy_comp <- re_combine %>% 
  filter(outcome=="amylsqrt") %>% 
  ggplot(aes(x=Beta_brain,y=Beta_serum,color=Class2))+
  geom_hline(yintercept = 0,linetype ="dashed",color="grey50")+
  geom_vline(xintercept = 0,linetype ="dashed",color="grey50")+
  geom_point()+
  geom_label_repel(aes(label = annotation),
                   size = 4,
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = "black",
                   max.overlaps = 30,
                   show.legend = FALSE)+
  scale_color_manual(values = class_col,na.value = "grey60",name="Class")+
  labs(x="Beta (Brain metabolites)",y="Beta (Serum metabolites)")+ #title = "Amyloid-beta"
  theme_classic()+
  theme(text = element_text(size = 18),
        legend.position = "none",
        plot.margin = margin(c(5.5,8,5.5,5.5),unit = "pt"))
print(p_amy_comp)
ggsave(p_amy_comp,
       filename = "figures/brain_serum/assoc_comparison_amyloid.pdf",
       width = 4,height = 4)

re_combine_tang <- re_combine %>% 
  filter(outcome=="tangsqrt")
head(re_combine_tang)
cor.test(re_combine_tang$Beta_brain,re_combine_tang$Beta_serum)


p_tang_comp <- re_combine %>% 
  filter(outcome=="tangsqrt") %>% 
  ggplot(aes(x=Beta_brain,y=Beta_serum,color=Class2))+
  geom_hline(yintercept = 0,linetype ="dashed",color="grey50")+
  geom_vline(xintercept = 0,linetype ="dashed",color="grey50")+
  geom_point()+
  geom_text_repel(aes(label = annotation),
                   size = 4,
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = "black",
                   max.overlaps = 20,
                   show.legend = FALSE)+
  scale_color_manual(values = class_col,na.value = "grey60",name="Class")+
  labs(x="Beta (Brain metabolites)",y="Beta (Serum metabolites)")+ #,title = "Tau tangles"
  theme_classic()+
  theme(text = element_text(size = 18),
        legend.position = "none",
        plot.margin = margin(c(5.5,8,5.5,5.5),unit = "pt"))
print(p_tang_comp)
ggsave(p_tang_comp,
       filename = "figures/brain_serum/assoc_comparison_tangle.pdf",
       width = 4,height = 4)


class_col = c("Acylcarnitines" = col_pal[1],
              "Aminoacids" = col_pal[2],
              "Biogenic Amines" = col_pal[3],
              "Glycerides" = col_pal[4],
              "Glycerophospholipids" = col_pal[5],
              "Sphingolipids" = col_pal[6],
              # "Cholesterol Esters" = col_pal[7],
              "Sugars" = col_pal[8])

p_leg <- ggplot(data.frame(Class = names(class_col), x = 1, y = 1),
                aes(x, y, color = Class)) +
  geom_point(size = 3) +
  scale_color_manual(values = class_col, name = "Class") +
  theme_void() +
  theme(legend.position = "right")+
  guides(color = guide_legend(nrow = 2))
p_leg

legend <- cowplot::get_legend(p_leg)
legend

pdf("figures/brain_serum/legend_comparison.pdf",height = 1.5,width = 6)
grid.draw(legend)
dev.off()

