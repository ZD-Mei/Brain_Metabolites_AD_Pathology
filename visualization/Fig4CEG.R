rm(list = ls())
library(tidyverse)
library(glmnet)
library(gtools)
library(openxlsx)
library(broom)
library(ggrepel)
library(vegan)
library(sjmisc)

select <- dplyr::select
int <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}

setwd("workpath")

# mets data -----
load("data/dat_combine_all.RData")


# snRNAseq data ------
snRNAseq <- read_rds("pastel_resources_data2share_snRNA_pseudoBulk_7majCellTypes.rds")
head(snRNAseq$ext[,1:10]);dim(snRNAseq$ext)

# filtered res -----
re_sn_met <- read.xlsx("snRNA/res_snRNA_mets_filter.xlsx")
head(re_sn_met);dim(re_sn_met)

(gene_list <- unique(re_sn_met$gene_id))

# ext ------
ext <- snRNAseq$ext %>% 
  rotate_df(.) %>% 
  select(any_of(gene_list)) %>% 
  rownames_to_column(var = "projid") %>% 
  mutate(projid=as.numeric(projid))
head(ext[,1:5]);dim(ext) #424/16835
summary(ext$ENSG00000238009)

## inverse normal transformation
ext[,2:ncol(ext)] <- map_df(ext[,2:ncol(ext)],int)

ext_gene_list <- colnames(ext)[2:ncol(ext)]
length(ext_gene_list) #16834

table(ext$projid %in% dat_all$projid) #315/109
dat_ext <- dat_all %>% 
  inner_join(ext,by="projid")
dim(dat_ext) #322
table(colSums(is.na(dat_ext)))

# inh ------
inh <- snRNAseq$inh %>% 
  rotate_df(.) %>% 
  select(any_of(gene_list)) %>% 
  rownames_to_column(var = "projid") %>% 
  mutate(projid=as.numeric(projid))
head(inh[,1:5]);dim(inh) #424/16468
summary(inh$ENSG00000241860)
## inverse normal transformation
inh[,2:ncol(inh)] <- map_df(inh[,2:ncol(inh)],int)

inh_gene_list <- colnames(inh)[2:ncol(inh)]
length(inh_gene_list) #16834

table(inh$projid %in% dat_all$projid) #102/322
dat_inh <- dat_all %>% 
  inner_join(inh,by="projid")
dim(dat_inh)


# oli ------
oli <- snRNAseq$oli %>% 
  rotate_df(.) %>% 
  select(any_of(gene_list)) %>% 
  rownames_to_column(var = "projid") %>% 
  mutate(projid=as.numeric(projid))
head(oli[,1:5]);dim(oli) #424/15860
summary(oli$ENSG00000241860)
## inverse normal transformation
oli[,2:ncol(oli)] <- map_df(oli[,2:ncol(oli)],int)

oli_gene_list <- colnames(oli)[2:ncol(oli)]
length(oli_gene_list) #16834

table(oli$projid %in% dat_all$projid) #102/322
dat_oli <- dat_all %>% 
  inner_join(oli,by="projid")
dim(dat_oli)


# ast ------
ast <- snRNAseq$ast %>% 
  rotate_df(.) %>% 
  select(any_of(gene_list)) %>% 
  rownames_to_column(var = "projid") %>% 
  mutate(projid=as.numeric(projid))
head(ast[,1:5]);dim(ast) # 424/16384
summary(ast$ENSG00000241860)
## inverse normal transformation
ast[,2:ncol(ast)] <- map_df(ast[,2:ncol(ast)],int)

ast_gene_list <- colnames(ast)[2:ncol(ast)]
length(ast_gene_list) #16834

table(ast$projid %in% dat_all$projid) #102/322
dat_ast <- dat_all %>% 
  inner_join(ast,by="projid")
dim(dat_ast)


# opc ------
opc <- snRNAseq$opc %>% 
  rotate_df(.) %>% 
  select(any_of(gene_list)) %>% 
  rownames_to_column(var = "projid") %>% 
  mutate(projid=as.numeric(projid))
head(opc[,1:5]);dim(opc) #424/15256
summary(opc$ENSG00000237491)
## inverse normal transformation
opc[,2:ncol(opc)] <- map_df(opc[,2:ncol(opc)],int)

opc_gene_list <- colnames(opc)[2:ncol(opc)]
length(opc_gene_list) #16834

table(opc$projid %in% dat_all$projid) #102/322
dat_opc <- dat_all %>% 
  inner_join(opc,by="projid")
dim(dat_opc)


# mic ------
mic <- snRNAseq$mic %>% 
  rotate_df(.) %>% 
  select(any_of(gene_list)) %>% 
  rownames_to_column(var = "projid") %>% 
  mutate(projid=as.numeric(projid))
head(mic[,1:5]);dim(mic) #424/14396
summary(mic$ENSG00000237491)
## inverse normal transformation
mic[,2:ncol(mic)] <- map_df(mic[,2:ncol(mic)],int)

mic_gene_list <- colnames(mic)[2:ncol(mic)]
length(mic_gene_list) #16834

table(mic$projid %in% dat_all$projid) #102/322
dat_mic <- dat_all %>% 
  inner_join(mic,by="projid")
dim(dat_mic)


# end ------
end <- snRNAseq$end %>% 
  rotate_df(.) %>% 
  select(any_of(gene_list)) %>% 
  rownames_to_column(var = "projid") %>% 
  mutate(projid=as.numeric(projid))
head(end[,1:5]);dim(end) #424/8413
summary(end$ENSG00000228794)
## inverse normal transformation
end[,2:ncol(end)] <- map_df(end[,2:ncol(end)],int)

end_gene_list <- colnames(end)[2:ncol(end)]
length(end_gene_list) #16834

table(end$projid %in% dat_all$projid) #102/322
dat_end <- dat_all %>% 
  inner_join(end,by="projid")
dim(dat_end)

# viz -----
head(re_sn_met)
unique(re_sn_met$Compounds)

i <- 1
for (i in seq_len(nrow(re_sn_met))) {
  celltype <- re_sn_met[i,"celltype"]
  gene_id <- re_sn_met[i,"gene_id"]
  com_id <- re_sn_met[i,"met"]
  compound <- re_sn_met[i,"Compounds"]
  gene_name <- re_sn_met[i,"gene_name"]
  data <- get(paste0("dat_",celltype))
  
  n_eff <- data %>%
    dplyr::filter(
      !is.na(.data[[com_id]]),
      !is.na(.data[[gene_id]])
    ) %>%
    nrow()
  
  
  p <- ggplot(data,aes(y=.data[[com_id]],x=.data[[gene_id]]))+
    geom_point(color="#87a2b3",size = 3,alpha=0.8)+
    geom_smooth(method = "lm")+
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0("n=",n_eff),
      hjust = -0.1, vjust = 1.1,
      size = 6
    ) +
    labs(x=paste0(gene_name," in ",celltype), y=compound)+
    theme_classic()+
    theme(text = element_text(size = 20))
  print(p)
  ggsave(p,filename = paste0("figures/snRNA_mets/scatter_",compound,"_",
                             celltype,"_",gene_name,".pdf"),
         width = 4,height = 3.5)
}



