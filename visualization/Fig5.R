rm(list = ls())
library(tidyverse)
library(data.table)
library(openxlsx)
library(sjmisc)
library(TwoSampleMR)
select <- dplyr::select

setwd("workpath")

load("data/dat_combine_all.RData")


# _OneSample significant ------
## mets -> amy ----
head(re_mets2amy_1sample)

re_mets2amy_1sample_sig <- 
  re_mets2amy_1sample %>% 
  filter(fdr_f < 0.1) %>% 
  filter(!duplicated(com_id,outcome)) %>%
  left_join(brain_targeted_info,by="com_id") %>% 
  arrange(Beta_f) %>% 
  mutate(
    fdr_pound=case_when(
      fdr_f < 0.1 ~ "#",
      TRUE ~ NA_character_),
    fdr_star=case_when(
      fdr_f >=0.1 & P_val_f <0.05 ~ "*",
      TRUE ~ NA_character_),
    sig.sign=case_when(
      fdr_f <0.1 ~ "#",
      fdr_f >=0.1 & P_val_f <0.05 ~ "*",
      P_val_f >=0.05 ~ "")
    )
re_mets2amy_1sample_sig
re_mets2amy_1sample_sig$Compounds <- 
  factor(re_mets2amy_1sample_sig$Compounds,
         levels = rev(c("PC(29:0)","LPC(20:3)","SM(30:1)","Glucose")))

p_mets2amy <- 
  ggplot(re_mets2amy_1sample_sig,aes(y=Compounds,x=Beta_f))+
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.6)+
  geom_errorbar(aes(xmin=Beta_f-1.96*SE_f,xmax=Beta_f+1.96*SE_f),
                linewidth=0.8,width=0)+
  geom_point(shape=22,fill="#3ca2e2",size=6,stroke=0.8)+
  geom_text(aes(label = Compounds),size=6,nudge_y = -0.3,na.rm=TRUE) +
  scale_y_discrete(position = "right") +
  coord_cartesian(
    xlim = c(-0.8, 0.8),
    clip = "off"
  ) +
  annotate(
    "segment",
    x = -0.85, xend = 0.85,
    y = 0.35, yend = 0.35,
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    linewidth = 0.8,
    color = "black"
  ) +
  theme_bw(base_size = 12) +
  theme(
    text=element_text(size=18,color="black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  )
p_mets2amy
ggsave("figures/mr/forest_mets_to_amy.pdf",p_mets2amy,
       width = 3,height = 3)


## mets -> tang -----
head(re_mets2tang_1sample)

re_mets2tang_1sample_sig <- 
  re_mets2tang_1sample %>% 
  filter(fdr_f < 0.1) %>% 
  filter(!duplicated(com_id,outcome)) %>% 
  left_join(brain_targeted_info,by="com_id") %>% 
  arrange(Beta_f) %>% 
  mutate(
    fdr_pound=case_when(
      fdr_f < 0.1 ~ "#",
      TRUE ~ NA_character_),
    fdr_star=case_when(
      fdr_f >=0.1 & P_val_f <0.05 ~ "*",
      TRUE ~ NA_character_),
    sig.sign=case_when(
      fdr_f <0.1 ~ "#",
      fdr_f >=0.1 & P_val_f <0.05 ~ "*",
      P_val_f >=0.05 ~ "")
  )
re_mets2tang_1sample_sig


p_mets2tang <-
  ggplot(re_mets2tang_1sample_sig,aes(y=Compounds,x=Beta_f))+
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.6)+
  geom_errorbar(aes(xmin=Beta_f-1.96*SE_f,xmax=Beta_f+1.96*SE_f),
                linewidth=0.8,width=0)+
  geom_point(shape=22,fill="#3ca2e2",size=6,stroke=0.8)+
  geom_text(aes(label = Compounds),size=6,nudge_y = -0.3,na.rm=TRUE) +
  scale_y_discrete(position = "right") +
  coord_cartesian(
    xlim = c(-0.8, 0.8),
    clip = "off"
  ) +
  annotate(
    "segment",
    x = -0.85, xend = 0.85,
    y = 0.35, yend = 0.35,
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    linewidth = 0.8,
    color = "black"
  ) +
  theme_bw(base_size = 12) +
  theme(
    text=element_text(size=18,color="black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  )
p_mets2tang
ggsave("figures/mr/heatmap_mets_to_tang.pdf",
       p_mets2tang,
       width = 3,height = 1.5)

## tang -> mets ----
head(re_tang2mets_1sample)
re_tang2mets_sig <- re_tang2mets_1sample %>% 
  filter(fdr_b < 0.1) %>% 
  arrange(Beta_b) %>% 
  mutate(
    fdr_pound=case_when(
      fdr_b < 0.1 ~ "#",
      TRUE ~ NA_character_),
    fdr_star=case_when(
      fdr_b >=0.1 & P_val_b <0.05 ~ "*",
      TRUE ~ NA_character_),
    sig.sign=case_when(
      fdr_b <0.1 ~ "#",
      fdr_b >=0.1 & P_val_b <0.05 ~ "*",
      P_val_b >=0.05 ~ "")
  )
re_tang2mets_sig

re_tang2mets_sig$Compounds <- 
  factor(re_tang2mets_sig$Compounds,
         levels = re_tang2mets_sig$Compounds)

p_tang2mets <- 
  ggplot(re_tang2mets_sig,aes(y=Compounds,x=Beta_b))+
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.6)+
  geom_errorbar(aes(xmin=Beta_b-1.96*SE_b,xmax=Beta_b+1.96*SE_b),
                linewidth=0.8,width=0)+
  geom_point(shape=22,fill="#3ca2e2",size=6,stroke=0.8)+
  geom_text(aes(label = Compounds),size=6,nudge_y = -0.3,na.rm=TRUE) +
  scale_y_discrete(position = "right") +
  coord_cartesian(
    xlim = c(-0.8, 0.8),
    clip = "off"
  ) +
  annotate(
    "segment",
    x = -0.85, xend = 0.85,
    y = 0.35, yend = 0.35,
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    linewidth = 0.8,
    color = "black"
  ) +
  theme_bw(base_size = 12) +
  theme(
    text=element_text(size=18,color="black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  )
p_tang2mets

ggsave("figures/mr/forest_tang_to_mets.pdf",p_tang2mets,
       width = 3,height = 1.5)


