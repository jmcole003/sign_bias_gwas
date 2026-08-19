#!/bin/bash

### Run binning and MAF threshold calculations


####################################################################################
#UKB

while read -r f; do

  for method in sig random; do

    #bins
    (
      cd /binned_out/UKB/bin

      Rscript /scripts/7-3_Bin_and_Calc_sign_bias.R \
        --file "$f" \
        --type MAF \
        --method "$method" \
        --mode bin

    ) > "/binned_out/UKB/logs/$(basename "$f").${method}.bin.log" 2>&1


    #thresholds
    for T in 0.001 0.01 0.1; do

      (
        cd /binned_out/UKB/threshold

        Rscript /scripts/7-3_Bin_and_Calc_sign_bias.R \
          --file "$f" \
          --type MAF \
          --method "$method" \
          --mode threshold \
          --low "$T"

      ) > "/binned_out/UKB/logs/$(basename "$f").${method}.thr_${T}.log" 2>&1

    done

  done

done < <(
  find /UKB/ash_out -type f \
    \( -name "*.txt" -o \
       -name "*.tsv" -o \
       -name "*.txt.gz" -o \
       -name "*.tsv.gz" \) | sort
)


####################################################################################
#AoU

while read -r f; do

  for method in sig random; do

    #bins
    (
      cd /binned_out/AOU/bin

      Rscript /scripts/7-3_Bin_and_Calc_sign_bias.R \
        --file "$f" \
        --type MAF \
        --method "$method" \
        --mode bin

    ) > "/binned_out/AOU/logs/$(basename "$f").${method}.bin.log" 2>&1


    #thresholds
    for T in 0.001 0.01 0.1; do

      (
        cd /binned_out/AOU/threshold

        Rscript /scripts/7-3_Bin_and_Calc_sign_bias.R \
          --file "$f" \
          --type MAF \
          --method "$method" \
          --mode threshold \
          --low "$T"

      ) > "/binned_out/AOU/logs/$(basename "$f").${method}.thr_${T}.log" 2>&1

    done

  done

done < <(
  find /AoU/ash_out -type f \
    \( -name "*.txt" -o \
       -name "*.tsv" -o \
       -name "*.txt.gz" -o \
       -name "*.tsv.gz" \) | sort
)


####################################################################################
#FinnGen
while read -r f; do

  for method in sig random; do

    #bins
    (
      cd /binned_out/FG/bin

      Rscript /scripts/7-3_Bin_and_Calc_sign_bias.R \
        --file "$f" \
        --type MAF \
        --method "$method" \
        --mode bin

    ) > "/binned_out/FG/logs/$(basename "$f").${method}.bin.log" 2>&1


    #thresholds
    for T in 0.001 0.01 0.1; do

      (
        cd /binned_out/FG/threshold

        Rscript /scripts/7-3_Bin_and_Calc_sign_bias.R \
          --file "$f" \
          --type MAF \
          --method "$method" \
          --mode threshold \
          --low "$T"

      ) > "/binned_out/FG/logs/$(basename "$f").${method}.thr_${T}.log" 2>&1

    done

  done

done < <(
  find /FG/ash_out -type f \
    \( -name "*.txt" -o \
       -name "*.tsv" -o \
       -name "*.txt.gz" -o \
       -name "*.tsv.gz" \) | sort
)