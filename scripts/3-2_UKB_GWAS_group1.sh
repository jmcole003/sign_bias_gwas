#!/bin/bash
#SBATCH -J ukb_gwas_group1                              
#SBATCH -N 1                 
#SBATCH -n 1                 
#SBATCH -c 16                  
#SBATCH -t 5:00:00           
#SBATCH --array=1-22            

# Chrom array index
CHR=${SLURM_ARRAY_TASK_ID}

# Run PLINK2 GWAS (group 1)
plink2 --bgen chr${CHR}.bgen ref-first \
  --sample chr${CHR}.sample \
  --keep /phenos/ukb_samples_group1.txt \
  --pheno /phenos/ukb_phenotypes.txt \
  --pheno-name AD T1D T2D SCZ ASTHMA height bmi weight mono_pct neutro_pct baso_pct wbc rbc mch \
  --covar /phenos/ukb.covariates.txt \
  --covar-name sex age age2 sex_age sex_age2 PC1-PC20 \
  --covar-variance-standardize \
  --maf 0.001 \
  --max-maf 0.999 \
  --vif 80 \
  --hwe 1e-10 midp keep-fewhet \
  --geno 0.05 \
  --glm hide-covar cols=+a1freq \
  --threads 16 \
  --memory 60000 \
  --out gwas_group1_chr${CHR}