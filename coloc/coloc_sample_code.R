
rm(list = ls())

library(tidyverse)
library(data.table)
library(openxlsx)
library(coloc)
select <- dplyr::select

setwd("workpath")

load("data/dat_combine_all.RData")
head(brain_targeted_info);dim(brain_targeted_info)
table(brain_targeted_info$Class)

# mQTL ------
re_mQTL_path <- "for_coloc/"
met_snp = "METSNP"
re_mqtl <- fread(paste0(re_mQTL_path,met_snp,".csv"))
head(re_mqtl);dim(re_mqtl) 

(chr <- str_extract(met_snp, "chr\\d{1,2}"))
(chr_number <- as.numeric(str_extract(chr, "\\d+")))

(met=str_extract(met_snp, "brain1_com\\d+"))
(met_name <- brain_targeted_info$Compounds[brain_targeted_info$com_id==met])

# eQTL ------
load("coloc/data_eqtl/rosmap_sneQTL_clean_rsID_match.RData")
head(re_eqtl_filtered);dim(re_eqtl_filtered) 
re_eqtl_filtered <- re_eqtl_filtered %>% 
  filter(!is.na(CHR)) %>% 
  filter(CHR==chr_number) 

celltype_gene_list <- unique(re_eqtl_filtered$celltype_gene)
length(celltype_gene_list) 
length(unique(re_eqtl_filtered$gene_symbol))

# coloc --------
re_coloc_path <- "coloc/mqtl_eqtl/"

no_overlap_list <- list()
for (cell_gene in celltype_gene_list) {
  re_eqtl <- re_eqtl_filtered %>% 
    filter(celltype_gene==cell_gene)
  head(re_eqtl);dim(re_eqtl)
  
  ## Check for SNP overlap 
  snp_list <- intersect(re_eqtl$ID, re_mqtl$ID)
  length(snp_list)
  
  if (length(snp_list) == 0) {
    message("No overlapping SNPs found, skipping coloc analysis.")
    no_overlap_list <- append(no_overlap_list, cell_gene)
  } else {
    print(cell_gene)
    merged_data <- re_mqtl %>% 
      inner_join(re_eqtl,by=c("SNP"="snps","ID","A1.ref","A2.ref","CHR","BP")) %>% 

      ## flip the beta if the alternative allele doesn't match
      mutate(
        flip = (ALLELE1 != ALT),  
        beta = ifelse(flip, -beta, beta),  # revise the direction of eQTL beta
        varbeta_mqtl = SE**2,
        varbeta_eqtl = se**2)
    
    ## mQTL data
    mqtl_coloc <- merged_data %>%
      select(SNP, CHR, GENPOS, BETA, varbeta_mqtl) 
    dim(mqtl_coloc)
    colnames(mqtl_coloc) <- c("snp", "chr", "position", "beta", "varbeta")
    mqtl_coloc_list <- as.list(mqtl_coloc)
    mqtl_coloc_list$type <- "quant"
    mqtl_coloc_list$sdY <- 1
    
    ## eQTL data
    eqtl_coloc <- merged_data %>%
      select(SNP, CHR, BP, beta, varbeta_eqtl)
    dim(eqtl_coloc)
    colnames(eqtl_coloc) <- c("snp", "chr", "position", "beta", "varbeta")
    eqtl_coloc_list <- as.list(eqtl_coloc)
    eqtl_coloc_list$type <- "quant"
    eqtl_coloc_list$sdY <- 1
    
    ## coloc 
    my.res <- coloc.abf(dataset1 = mqtl_coloc_list, dataset2 = eqtl_coloc_list)
    max_pp_snp = my.res$results[which.max(my.res$results$SNP.PP.H4),]
    
    re_coloc <- data.table(
      met_snp = met_snp,
      Compounds = met_name,
      celltype_gene = cell_gene,
      PP.H0.abf = my.res$summary["PP.H0.abf"],
      PP.H1.abf = my.res$summary["PP.H1.abf"],
      PP.H2.abf = my.res$summary["PP.H2.abf"],
      PP.H3.abf = my.res$summary["PP.H3.abf"],
      PP.H4.abf = my.res$summary["PP.H4.abf"],
      coloc_region = paste0(chr_number,":",min(merged_data$GENPOS),
                            "-",max(merged_data$GENPOS)),
      nsnps_coloc = my.res$summary["nsnps"],
      snp_max_pp = max_pp_snp$snp,
      snp_max_pp_chr = chr_number,
      snp_max_pp_position = max_pp_snp$position,
      max_pp_val = max_pp_snp$SNP.PP.H4,
      
      
      mQTL_nsnps_raw = nrow(re_mqtl),
      mQTL_p_causal = merged_data$PVALUE[which(merged_data$SNP==max_pp_snp$snp)],

      # mQTL SNP with minimum p value of all
      mQTL_min_p_raw = min(re_mqtl$PVALUE,na.rm = T),
      mQTL_min_p_raw_snp = re_mqtl$SNP[which.min(re_mqtl$PVALUE)],

      # mQTL SNP with minimum p value for coloc 
      mQTL_min_p_coloc = min(merged_data$PVALUE),
      mQTL_min_p_coloc_snp = merged_data$SNP[which.min(merged_data$PVALUE)],
      
      eQTL_nsnps_raw = nrow(re_eqtl),
      eQTL_p_causal = merged_data$pvalue[which(merged_data$SNP==max_pp_snp$snp)],

      # mQTL SNP with minimum p value of all
      eQTL_min_p_raw = min(re_eqtl$pvalue),
      eQTL_min_p_raw_snp = re_eqtl$snps[which.min(re_eqtl$pvalue)],
      
      # mQTL SNP with minimum p value for coloc 
      eQTL_min_p_coloc = min(merged_data$pvalue),
      eQTL_min_p_coloc_snp = merged_data$SNP[which.min(merged_data$pvalue)])
      
      # Save results
      fwrite(re_coloc, paste0(re_coloc_path,"coloc_",
                            met_snp, "_", cell_gene, ".csv"))
  }
}

if (length(no_overlap_list) > 0) {
  no_overlap_file <- 
    paste0("coloc/mqtl_eqtl_nosnps/no_overlap_summary_",met_snp,".txt")
  write_lines(no_overlap_list, no_overlap_file)
}


