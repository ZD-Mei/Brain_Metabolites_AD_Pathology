import subprocess

# File paths
phenotype_list_file = "phenotype_list_brain_mets.txt"
log_dir = "logDir/"
out_dir = "mod_out/"
pred_file = "step1_output_mets_pred.list"
covar_file = "covariates_brain_mets_mod1.txt"
pheno_file = "phenotype_brain_mets.txt"

# Read phenotype_list.txt
with open(phenotype_list_file, "r") as f:
    phenotypes = [line.strip() for line in f if line.strip()]

# Iterate through each phenotype and submit the job
for pheno in phenotypes:
    print(f"Submitting job for phenotype: {pheno}")

    sbatch_command = [
        "sbatch",
        "--job-name", f"regenie_step2_{pheno}",
        "--export=ALL",
        "--cpus-per-task=8",
        "--time=0-12:00",
        "--partition=sapphire",
        "--mem=500000",
        f"--output={log_dir}regenie_step2_{pheno}_%j.out",
        f"--error={log_dir}regenie_step2_{pheno}_%j.err",
        "--wrap",
        f"""
        module load Anaconda/5.0.1-fasrc01
        source activate regenie_env
        regenie \
          --step 2 \
          --bed /data_path/combined_qc \
          --covarFile {covar_file} \
          --covarColList age_death,msex,pmi,PC1,PC2,PC3,PC4 \
          --phenoFile {pheno_file} \
          --phenoCol {pheno} \
          --bsize 1000 \
          --threads 8 \
          --gz \
          --pred {pred_file} \
          --out {out_dir}assoc_{pheno}
        """
    ]

    subprocess.run(sbatch_command, check=True)
