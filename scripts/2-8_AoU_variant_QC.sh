#!/bin/bash
## QC AoU variants

pgen_dir="./plink"
keep_file="./wb_samples_filtered.txt"
rsid_file="./rsids_from_ukb_list.txt"

#QC chr 1-22
for chr in $(seq 1 22); do

  plink2 \
    --pfile "${pgen_dir}/chr${chr}/chr${chr}" \
    --keep "$keep_file" \
    --hwe 1e-12 midp keep-fewhet \
    --geno 0.05 \
    --extract "$rsid_file" \
    --make-bed \
    --out "${pgen_dir}/chr${chr}/chr${chr}.filtered" \
    --allow-extra-chr

done