rm(list = ls())
library(tidyverse)
library(data.table)
library(openxlsx)


setwd("workpath")

select <- dplyr::select


load("data/dat_combine_all.RData")

dat_all$ser_visit
summary(dat_all$age_at_visit)
summary(dat_all$age_death)

dat_all$age_diff <- dat_all$age_death - dat_all$age_at_visit
summary(dat_all$age_diff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.558   5.002   7.110   7.637   9.822  16.564

dat_brain_serum <- dat_all %>% 
  filter(!is.na(brain1_com1)&!is.na(serum1_com1))
dim(dat_brain_serum)
colSums(is.na(dat_brain_serum))

# age difference -------
ggplot(dat_all,aes(x=age_diff))+
  geom_histogram(
    bins = 30,                
    color = "white",            
    fill = "#69b3a2",           
    alpha = 0.8                 
  ) +
  labs(x = "Years", y = "Count") + 
  theme_classic() +           
  theme(
    text = element_text(size = 18),
    plot.title = element_text(hjust = 0.5) # 居中标题
  )
ggsave("figures/brain_serum/fig6a_age_difference_histogram.pdf",
       width = 4,height = 3)


# density plot -------
line_positions <- dat_age %>%
  group_by(age_group) %>%
  summarise(vline = median(age))
p_age_density <- 
  ggplot(dat_age,aes(x=age,fill=age_group))+
  geom_vline(data = line_positions, aes(xintercept = vline, color = age_group),
             linetype = "dashed", size = 1) +
  geom_density(
    color = "white",            
    alpha = 0.6  
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e"),
                    labels = c("Age at blood sample collection",
                               "Age at death (brain sample collection)"),
                    name = "") +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e"),guide = "none")+
  labs(x = "Years", y = "Density") +
  theme_classic() +       
  theme(
    text = element_text(size = 18),
    plot.title = element_text(hjust = 0.5), # 居中标题
    legend.position = "none",
    plot.margin = unit(c(5.5,9,5.5,5.5),"pt")
  )
p_age_density
ggsave(p_age_density,
       filename = "figures/brain_serum/fig6a_age_density_plot.pdf",
       width = 4,height = 3)

p_age_density_lgd <- p_age_density+
  theme(legend.position = "right",
        legend.text  = element_text(size = 18))
p_age_density_legend <- ggpubr::get_legend(p_age_density_lgd)

pdf("figures/brain_serum/fig6a_age_density_plot_legend.pdf",
    height = 1,width = 5)
grid.draw(p_age_density_legend)
dev.off()

