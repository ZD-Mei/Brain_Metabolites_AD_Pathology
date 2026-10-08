rm(list = ls())
library(tidyverse)
library(data.table)
library(openxlsx)
library(CMAverse)
# library(glue)

setwd("workpath")

# load data ------
load("data/dat_combine_all.RData")
dat_brain <- dat_all %>% 
  filter(!is.na(brain1_com1))

# sel-define function -------
run_cmest_analysis <- 
  function(com_id, exposure, mediator, outcome,  
           yreg = "linear") {
  
    covs <- c("age_death", "msex", "educ", "pmi")
    dat_mediation <- dat_brain %>% 
      select(com_id,exposure,mediator,outcome,covs) %>% 
      tidyr::drop_na()
    a0 <- mean(dat_mediation[[exposure]], na.rm = TRUE)
    a1 <- a0 + sd(dat_mediation[[exposure]], na.rm = TRUE)
    mvalue <- mean(dat_mediation[[mediator]],na.rm = TRUE)
    
    tryCatch({
      set.seed(123)
      est <- cmest(
        data = dat_mediation, model = "rb",
        exposure = exposure, mediator = mediator, outcome = outcome,
        basec = covs, EMint = FALSE,
        mreg = list("linear"), yreg = yreg,
        multimp = FALSE, astar = a0, a = a1,
        mval = list(mvalue), #yval = yval,
        estimation = "paramfunc", inference = "delta"
      )
      
      results <- tibble(
        exposure = exposure,
        mediator = mediator,
        outcome = outcome,
        n = nrow(dat_mediation),
        effect = names(est$effect.pe),
        Estimate = est$effect.pe,
        SE = est$effect.se,
        CI.low = est$effect.ci.low,
        CI.high = est$effect.ci.high,
        P.val = est$effect.pval
      )
    }, error = function(e) {
      message("Error for ", com_id, ": ", conditionMessage(e))
      tibble(
        exposure = exposure,
        mediator = mediator,
        outcome = outcome,
        n = nrow(dat_mediation),
        effect = NA_character_,
        Estimate = NA_real_,
        SE = NA_real_,
        CI.low = NA_real_,
        CI.high = NA_real_,
        P.val = NA_real_
      )
    })
  }


com_id = "METID"

analysis_configs <- tibble(
  exposure = c(com_id, com_id, com_id, 
               "amylsqrt", "amylsqrt", "tangsqrt"),
  mediator = c("amylsqrt", "amylsqrt", "tangsqrt",
               com_id, com_id, com_id),
  outcome = c("tangsqrt", "cogng_demog_slope", "cogng_demog_slope", 
              "tangsqrt", "cogng_demog_slope", "cogng_demog_slope"),
  yreg = c("linear", "linear", "linear", "linear", "linear", "linear")
)


# Analyze and combine results 
final_results <- pmap_dfr(
  analysis_configs, 
  ~run_cmest_analysis(com_id,..1, ..2, ..3, ..4)
)
dim(final_results)

# save final results
fwrite(final_results, 
       paste0("results/mediation/med_framework_cmaverse/re_med_framework_",
              com_id,".csv"))

