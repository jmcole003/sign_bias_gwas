#!/bin/bash

## AoU GWAS
pheno_dir="phenos"
geno_dir="plink"
out_dir="gwas_out"

covar="${pheno_dir}/aou_covariates.centered.tsv"
covar_plink="${pheno_dir}/aou_covariates.centered.plink.tsv"

threads=8
max_jobs=5

####################################################################################
## covariates

#add FID/IID for plink

awk -F'\t' -v OFS='\t' '
NR==1{
  printf "FID\tIID"
  for(i=2;i<=NF;i++) printf "\t%s",$i
  printf "\n"
  next
}
{
  printf "%s\t%s",$1,$1
  for(i=2;i<=NF;i++) printf "\t%s",$i
  printf "\n"
}
' "$covar" > "$covar_plink"


####################################################################################
## binary phenotypes

#recode 
binary_traits=("asthma" "t1d" "t2d" "schizophrenia" "alzheimers")

for trait in "${binary_traits[@]}"; do

  infile="${pheno_dir}/aou_${trait}.plink.tsv"
  outfile="${pheno_dir}/aou_${trait}.linear.plink.tsv"

  if awk -F'\t' 'NR>1 && $3==0 {found=1} END{exit(found ? 0 : 1)}' "$infile"; then

    awk -F'\t' -v OFS='\t' '
    NR==1{print; next}
    {
      if($3==0) $3=3
      if($3==1) $3=4
      print
    }
    ' "$infile" > "$outfile"

  else

    awk -F'\t' -v OFS='\t' '
    NR==1{print; next}
    {
      if($3==1) $3=3
      if($3==2) $3=4
      print
    }
    ' "$infile" > "$outfile"

  fi

done


####################################################################################
## GWAS

traits=("height"
        "weight"
        "bmi"
        "monocyte_percentage"
        "basophil_percentage"
        "neutrophil_percentage"
        "white_blood_cell_count"
        "red_blood_cell_count"
        "mean_corpuscular_hemoglobin"
        "asthma"
        "t1d"
        "t2d"
        "schizophrenia"
        "alzheimers")


run_gwas () {

  trait=$1

  mkdir -p "${out_dir}/${trait}"

  if [[ " ${binary_traits[*]} " == *" ${trait} "* ]]; then
    pheno="${pheno_dir}/aou_${trait}.linear.plink.tsv"
  else
    pheno="${pheno_dir}/aou_${trait}.plink.tsv"
  fi


  #chr 1-22
  for chr in $(seq 1 22); do

    plink2 \
      --bfile "${geno_dir}/chr${chr}/chr${chr}.filtered" \
      --pheno "$pheno" \
      --covar "$covar_plink" \
      --covar-variance-standardize \
      --glm hide-covar omit-ref \
      --vif 80 \
      --threads "$threads" \
      --out "${out_dir}/${trait}/chr${chr}" \
      > "${out_dir}/logs/${trait}.chr${chr}.log" 2>&1

  done


  #combine 
  case "$trait" in
    height) outname="Height_AoU" ;;
    weight) outname="Weight_AoU" ;;
    bmi) outname="BMI_AoU" ;;
    monocyte_percentage) outname="Monocyte_percentage_AoU" ;;
    basophil_percentage) outname="Basophil_percentage_AoU" ;;
    neutrophil_percentage) outname="Neutrophil_percentage_AoU" ;;
    white_blood_cell_count) outname="White_blood_cell_count_AoU" ;;
    red_blood_cell_count) outname="Red_blood_cell_count_AoU" ;;
    mean_corpuscular_hemoglobin) outname="Mean_corpuscular_hemoglobin_AoU" ;;
    asthma) outname="Asthma_AoU" ;;
    t1d) outname="Type_1_diabetes_AoU" ;;
    t2d) outname="Type_2_diabetes_AoU" ;;
    schizophrenia) outname="Schizophrenia_AoU" ;;
    alzheimers) outname="Alzheimers_AoU" ;;
  esac


  awk 'FNR==1 && NR!=1 {next} {print}' \
    "${out_dir}/${trait}"/chr*.glm.linear \
    > "${out_dir}/${outname}.tsv"

  gzip -f "${out_dir}/${outname}.tsv"
}


####################################################################################
## run

for trait in "${traits[@]}"; do

  while [ "$(jobs -r | wc -l)" -ge "$max_jobs" ]; do
    sleep 2
  done

  run_gwas "$trait" &

done