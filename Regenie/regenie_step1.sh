
module load Anaconda/5.0.1-fasrc01
source activate regenie_env


regenie \
  --step 1 \
  --bed genotyped_for_mets_step1 \
  --covarFile covariates_brain_mets.txt \
  --phenoFile phenotype_brain_mets.txt \
  --bsize 1000 \
  --threads 8 \
  --lowmem \
  --out step1_output_mets
  
