#!/bin/bash

## Get AoU genotype data (ACAF)
CDR_STORAGE_PATH="${CDR_STORAGE_PATH:-gs://fc-aou-datasets-controlled/v8}"

src="${CDR_STORAGE_PATH}/wgs/short_read/snpindel/acaf_threshold/pgen"
vat="${CDR_STORAGE_PATH}/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz"

out="${PWD}/plink"
mkdir -p "$out"

#genotypes
for chr in $(seq 1 22); do
  mkdir -p "${out}/chr${chr}"
  
  gsutil -m -u "$GOOGLE_PROJECT" cp -r \
    "${src}/chr${chr}/"* \
    "${out}/chr${chr}/"
done


#VAT rsids
gsutil -u "$GOOGLE_PROJECT" cat "$vat" | \
  gzip -dc | \
  awk -F'\t' '
    BEGIN{OFS="\t"}
    
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="contig") c=i
        if($i=="position") p=i
        if($i=="dbsnp_rsid") r=i
      }
      print "contig","position","dbsnp_rsid"
      next
    }
    
    $r ~ /^rs[0-9]+$/ {
      print $c,$p,$r
    }
  ' > vat_rsids.clean.tsv


echo "Done"