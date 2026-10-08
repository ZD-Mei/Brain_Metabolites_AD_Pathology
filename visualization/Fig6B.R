rm(list = ls())

library(data.table)
library(tidyverse)
library(openxlsx)
library(reshape2)
library(grid)
library(gridExtra)
library(missRanger)
library(scales)


setwd("workpath")

dir <- "data_dir"
brain_fname <- paste0(dir, "Biocrates.Tissue.data_20240315.xlsx")
serum_fname <- paste0(dir, "Biocrates.Serum.data_20240315.xlsx")

brain <- as.data.frame(read.xlsx(brain_fname))
serum <- as.data.frame(read.xlsx(serum_fname))
head(colnames(brain),20);dim(brain)
head(colnames(serum),20);dim(serum) 

table(brain$Class)

brain_targeted_info_all <- brain[,1:16]
write.xlsx(brain_targeted_info_all,"data/brain_targeted_info_all.xlsx")
serum_targeted_info_all <- serum[,1:16]
write.xlsx(serum_targeted_info_all,"data/serum_targeted_info_all.xlsx")


# Missing pattern before NA removal -----
## brain ------
dat_brain <- brain[,c(1,17:ncol(brain))] %>% 
  pivot_longer("46291609":"50409532",
               names_to = "Samples",
               values_to = "metabolite_abundance")
head(dat_brain);dim(dat_brain)
length(unique(dat_brain$Compounds)) #409
length(unique(dat_brain$Samples)) #458

brain_detect_rate <- data.table(
  Compounds = brain$Compounds,
  detect_rate_brain = 100 - rowMeans(is.na(brain[,c(17:ncol(brain))]))*100
)
head(brain_detect_rate)

## serum ------
colnames(serum)
dat_serum <- serum[,c(1,17:ncol(serum))] %>% 
  pivot_longer("78905055":"17260313",
               names_to = "Samples",
               values_to = "metabolite_abundance")
head(dat_serum);dim(dat_serum)
length(unique(dat_serum$Compounds)) #409
length(unique(dat_serum$Samples)) #488

table(colnames(dat_brain) %in% colnames(dat_serum))

serum_detect_rate <- data.table(
  Compounds = serum$Compounds,
  detect_rate_serum = 100 - rowMeans(is.na(serum[,c(17:ncol(serum))]))*100
)
head(serum_detect_rate)

## missing rate comparison (brain vs. serum) ------
table(brain_detect_rate$Compounds %in% serum_detect_rate$Compounds)
table(serum_detect_rate$Compounds %in% brain_detect_rate$Compounds)

detect_rate <- brain_detect_rate %>% 
  full_join(serum_detect_rate,by="Compounds") %>% 
  mutate(detection_group=case_when(
    detect_rate_brain>70 & detect_rate_serum>70 ~ "both",
    detect_rate_brain>70 & detect_rate_serum<=70 ~ "brain",
    detect_rate_brain<=70 & detect_rate_serum>70 ~ "serum",
    detect_rate_brain<=70 & detect_rate_serum<=70 ~ "neither",
    TRUE ~ NA_character_
  ))
head(detect_rate)
detect_rate$detection_group <- 
  factor(detect_rate$detection_group,
         levels = c("both","serum","brain","neither"))
table(detect_rate$detection_group)

nejm_colors <- ggsci::pal_nejm()(8)
show_col(nejm_colors)

p_detect_rate <- 
  ggplot(detect_rate,aes(x=detect_rate_brain,y=detect_rate_serum,
                         color=detection_group))+
  geom_vline(xintercept = 70,linetype = "dashed")+
  geom_hline(yintercept = 70,linetype = "dashed")+
  scale_color_manual(
    values = c("#bc3c29","#0072b5","#20854e","#e18727"),
    labels = c("Both brain and serum (n=156)",
               "Serum only (n=114)",
               "Brain only (n=21)",
               "Neither brain nor serum (n=118)")
  )+
  geom_point()+
  scale_x_continuous(position = "top")+
  scale_y_reverse()+
  labs(x="Detection rate \nin brain samples (%)",
       y="Detection rate \nin serum samples (%)")+
  theme_bw() +
  theme(
    text = element_text(size = 18),
    panel.grid = element_blank(),
    legend.position = "none")
p_detect_rate

ggsave(filename = "figures/fig6b_detection_rate_brain_serum_before.pdf",
       p_detect_rate,
       width = 4,height = 4)

p_detect_rate_lgd <- p_detect_rate+
  theme(legend.position = "bottom",
        legend.title.align = 0.5
  )+
  guides(color = guide_legend(nrow = 2,
                              title.position = "top",
                              title = "Detected metabolites",
                              override.aes = list(size=3)))
p_detect_rate_lgd
p_detect_rate_legend <- ggpubr::get_legend(p_detect_rate_lgd)

pdf("figures/fig6b_detection_rate_brain_serum_before_legend.pdf",
    width = 8,height = 4)
grid.draw(p_detect_rate_legend)
dev.off()

