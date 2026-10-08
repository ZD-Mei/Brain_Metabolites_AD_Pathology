rm(list = ls())
library(tidyverse)
library(data.table)
library(openxlsx)
library(readr)
library(sjmisc)
library(RColorBrewer)
library(ggsci)
library(scales)
library(ggrepel)
library(qqman)

setwd("workpath")


load("data/dat_combine_all.RData")
head(brain_targeted_info)

# Mets list ------
mets_list <- read_lines("lists/ad_path_mets_list.txt")
head(mets_list);length(mets_list)

brain_targeted_info
mets_sel <- c("SM(34:1)","PC(36:3)","Asparagine","Tyrosine","Phenylalanine")
mets_sel_id <- brain_targeted_info %>% 
  filter(Compounds %in% mets_sel) %>% 
  pull(com_id)
mets_sel_id

# extract locus  ----
(met <- mets_sel_id[1])
locus_list <- list()
index = 1
for (met in mets_sel_id){
  met_name <- brain_targeted_info %>% 
    filter(com_id==met) %>% pull(Compounds)
  independent_snps <- 
    read.xlsx(paste0("gwas/brainMetabGWAS/resClump/gwas_",met,"_",
                     met_name,"_clump_filtered.xlsx"))
  
  re_gwas <- fread(paste0("gwas/brainMetabGWAS/reGWAS_filtered/assoc_",
                          met,".filtered.regenie.gz"))
  head(re_gwas);dim(re_gwas) # 5060574
  re_gwas$PVALUE <- 10^(-re_gwas$LOG10P)
  
  # extract peak area -----
  independent_snps
  window_size <- 500000 # 500kb
  # i=1
  for (i in seq_len(nrow(independent_snps))) {
    target_snp = independent_snps$SNP[i]
    target_chr = independent_snps$chr.exposure[i]
    target_pos = independent_snps$pos.exposure[i]
    gwas_subset <- re_gwas %>%
      filter(CHROM == target_chr 
             & GENPOS >= (target_pos - window_size) 
             & GENPOS <= (target_pos + window_size)) %>%
      mutate(Met = met_name) %>% 
      arrange(GENPOS)
    fwrite(gwas_subset,
           paste0("gwas/brainMetabGWAS/combined_mhp/snp_list/",
                  met,"_",met_name,"_chr",target_chr,"_",target_pos,
                  "_",target_snp,".csv"))
    locus_list[[index]] <- gwas_subset
    index <- index + 1
  }
}

re_all_locus <- bind_rows(locus_list)
head(re_all_locus);dim(re_all_locus) #67565
table(duplicated(re_all_locus$ID))
# re_dup <- re_all_locus %>% 
#   filter(duplicated(ID))
re_all_locus_dedup <- re_all_locus %>% 
  arrange(desc(LOG10P)) %>% 
  filter(!duplicated(ID)) %>% 
  arrange(Met,CHROM,GENPOS)
head(re_all_locus_dedup);dim(re_all_locus_dedup) #64452

re_gwas_null <- 
  fread("gwas/brainMetabGWAS/reGWAS_filtered/assoc_com1.filtered.regenie.gz") %>% 
  filter(!ID %in% re_all_locus_dedup$ID) %>% 
  mutate(LOG10P = ifelse(LOG10P>=5,4,LOG10P),
         PVALUE = 10^(-LOG10P))
head(re_gwas_null);dim(re_gwas_null)

dat_comb_mhplot <- 
  bind_rows(re_gwas_null,re_all_locus_dedup) %>% 
  arrange(CHROM,GENPOS)
head(dat_comb_mhplot);dim(dat_comb_mhplot)

# mh plot ----
## manhattan function definition
mh_plot_comb <- 
  function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
            col = c("gray80","gray60"), chrlabs = NULL, 
            suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), 
            highlight1 = NULL, highlight2 = NULL,highlight3 = NULL,
            highlight4 = NULL, highlight5 = NULL,highlight6 = NULL,
            logp = TRUE, ...) {
    CHR = BP = P = index = NULL
    if (!(chr %in% names(x)))
      stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x)))
      stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x)))
      stop(paste("Column", p, "not found!"))
    if (!(snp %in% names(x)))
      warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!is.numeric(x[[chr]]))
      stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]]))
      stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]]))
      stop(paste(p, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
    if (!is.null(x[[snp]]))
      d = transform(d, SNP = x[[snp]])
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP), ]
    if (logp) {
      d$logp <- -log10(d$P)
    }
    else {
      d$logp <- d$P
    }
    d$pos = NA
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
      ind = ind + 1
      d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
      options(scipen = 999)
      d$pos = d$BP/1e+06
      ticks = floor(length(d$pos))/2 + 1
      xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
      labs = ticks
    }
    else {
      lastbase = 0
      ticks = NULL
      for (i in unique(d$index)) {
        if (i == 1) {
          d[d$index == i, ]$pos = d[d$index == i, ]$BP
        }
        else {
          lastbase = lastbase + tail(subset(d, index ==
                                              i - 1)$BP, 1)
          d[d$index == i, ]$pos = d[d$index == i, ]$BP +
            lastbase
        }
        ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + 
                            max(d[d$CHR == i, ]$pos))/2 + 1)
      }
      xlabel = "Chromosome"
      labs <- unique(d$CHR)
    }
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
                     las = 1, pch = 20, xlim = c(xmin, xmax), 
                     ylim = c(0,ceiling(max(d$logp))), xlab = xlabel, 
                     ylab = expression(-log10))
    dotargs <- list(...)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
                                              names(dotargs)]))
    if (!is.null(chrlabs)) {
      if (is.character(chrlabs)) {
        if (length(chrlabs) == length(labs)) {
          labs <- chrlabs
        }
        else {
          warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
        }
      }
      else {
        warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
      }
    }
    if (nchr == 1) {
      axis(1, ...)
    }
    else {
      axis(1, at = ticks, labels = labs, ...)
    }
    col = rep(col, max(d$CHR))
    if (nchr == 1) {
      with(d, points(pos, logp, pch = 20, col = col[1], ...))
    }
    else {
      icol = 1
      for (i in unique(d$index)) {
        with(d[d$index == unique(d$index)[i], ], 
             points(pos, logp, col = col[icol], pch = 20, ...))
        icol = icol + 1
      }
    }
    if (!is.null(suggestiveline)) {
      abline(h = suggestiveline, col = "grey",lty = "dashed")
    }
    if (!is.null(genomewideline)) {
      abline(h = genomewideline, col = "red")
    }
    if (!is.null(highlight1)) {
      if (any(!(highlight1 %in% d$SNP)))
        warning("You're trying to highlight1 SNPs that don't exist in your results.")
      d.highlight1 = d[which(d$SNP %in% highlight1), ]
      with(d.highlight1, points(pos, logp, col = "#BC3C29CC", pch = 20,
                                ...))
    }
    if (!is.null(highlight2)) {
      if (any(!(highlight2 %in% d$SNP)))
        warning("You're trying to highlight2 SNPs that don't exist in your results.")
      d.highlight2 = d[which(d$SNP %in% highlight2), ]
      with(d.highlight2, points(pos, logp, col = "#0072B5CC", pch = 20,
                                ...))
    }
    if (!is.null(highlight3)) {
      if (any(!(highlight3 %in% d$SNP)))
        warning("You're trying to highlight3 SNPs that don't exist in your results.")
      d.highlight3 = d[which(d$SNP %in% highlight3), ]
      with(d.highlight3, points(pos, logp, col = "#E18727CC", pch = 20,
                                ...))
    }
    if (!is.null(highlight4)) {
      if (any(!(highlight4 %in% d$SNP)))
        warning("You're trying to highlight4 SNPs that don't exist in your results.")
      d.highlight4 = d[which(d$SNP %in% highlight4), ]
      with(d.highlight4, points(pos, logp, col = "#7876B1CC", pch = 20,
                                ...))
    }
    if (!is.null(highlight5)) {
      if (any(!(highlight5 %in% d$SNP)))
        warning("You're trying to highlight5 SNPs that don't exist in your results.")
      d.highlight5 = d[which(d$SNP %in% highlight5), ]
      with(d.highlight5, points(pos, logp, col = "#6F99ADCC", pch = 20,
                                ...))
    }
    if (!is.null(highlight6)) {
      if (any(!(highlight6 %in% d$SNP)))
        warning("You're trying to highlight6 SNPs that don't exist in your results.")
      d.highlight6 = d[which(d$SNP %in% highlight6), ]
      with(d.highlight6, points(pos, logp, col = "#EE4C97CC", pch = 20,
                                ...))
    }
  }

(col_pal=pal_nejm("default",alpha=0.8)(8))
scales::show_col(col_pal)

table(re_all_locus_dedup$Met)
asn_set <- re_all_locus_dedup %>% filter(Met=="Asparagine") %>% pull(ID)
tyr_set <- re_all_locus_dedup %>% filter(Met=="Tyrosine") %>% pull(ID)
phy_set <- re_all_locus_dedup %>% filter(Met=="Phenylalanine") %>% pull(ID)
pc_36_3_set <- re_all_locus_dedup %>% filter(Met=="PC(36:3)") %>% pull(ID)
sm_34_1_set <- re_all_locus_dedup %>% filter(Met=="SM(34:1)") %>% pull(ID)

head(dat_comb_mhplot)

png("figures/combined_mh_plot/combined_mhplot.png",
    height = 10,width = 20,units = "cm",res = 300)
mh_plot_comb(dat_comb_mhplot, 
             chr="CHROM",bp="GENPOS",p = "PVALUE",snp = "ID", 
             genomewideline = NULL,cex = 0.6,
             highlight1 = asn_set,
             highlight2 = tyr_set,
             highlight3 = phy_set,
             highlight4 = pc_36_3_set,
             highlight5 = sm_34_1_set
)
dev.off()



