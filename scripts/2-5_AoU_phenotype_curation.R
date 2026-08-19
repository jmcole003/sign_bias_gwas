## AoU phenotype/sample QC

library(data.table)

bucket <- Sys.getenv("WORKSPACE_BUCKET")

#calculate moments
calc_moments <- function(x){
  
  x <- x[!is.na(x)]
  mx <- mean(x)
  mu2 <- mean((x-mx)^2)
  mu3 <- mean((x-mx)^3)
  
  skew <- mu3/(mu2^(3/2))
  
  out <- c(mu0=1,
           mu1=0,
           mu2=mu2,
           mu3=mu3,
           mu3_std=skew)
  
  return(out)
}


#write phenotype files
write_pheno <- function(dat, pheno, ids, outfile){
  
  tmp <- data.frame(person_id=trimws(as.character(dat$person_id)),
                    value=dat[[pheno]])
  tmp <- tmp[tmp$person_id %in% ids,]
  tmp <- tmp[!is.na(tmp$value),]
  
  write.table(tmp, outfile,
              sep="\t", row.names=FALSE, quote=FALSE)
  
  #plink
  pl <- data.frame(FID=tmp$person_id,
                   IID=tmp$person_id,
                   PHENO=tmp$value)
  pfile <- sub(".tsv", ".plink.tsv", outfile, fixed=TRUE)
  write.table(pl, pfile,
              sep="\t", row.names=FALSE, quote=FALSE)
  return(tmp)
}


#center covariates
center_covs<- function(dat, covars){
  x <- dat
  for(v in covars){
    mn <- mean(x[[v]], na.rm=TRUE)
    x[[v]] <- x[[v]]-mn
  }
  return(x)
}


####################################################################################
## download inputs

files <- c(
  "aou_gwas/phecodeX_files/mcc2_phecodex_table_v8.csv",
  "aou_gwas/ancestry_table.csv",
  "aou_gwas/demographics_table.csv",
  "aou_gwas/participant_PCs.csv",
  "aou_gwas/samples_relatedness_flagged_samples.tsv",
  "aou_gwas/aou.v8.classified-rf.wb-classify.report-any.fam",
  "aou_gwas/pheno/lab_measures_all_v2.tsv",
  "aou_gwas/pheno/physical_measures_hwb.tsv"
)

files <- paste0(bucket, "/", files)

system2("gsutil", c("-m", "cp", files, "."))


####################################################################################
## Read data

#binary phenotypes
bcols <- c(
  "person_id",
  "sex",
  "RE_475",
  "EM_202.1",
  "EM_202.2",
  "MB_287.1",
  "NS_328.11"
)

binary <- fread(
  "mcc2_phecodex_table_v8.csv",
  select=bcols
)

binary[, person_id := trimws(as.character(person_id))]

fwrite(
  binary,
  file.path(phdir, "binary_phenos.csv")
)

#demographics
demo <- fread(
  "demographics_table.csv",
  select=1:3
)

setnames(demo, c("person_id","age","sex"))
demo[, person_id := trimws(as.character(person_id))]


#PCs
pcs <- fread(
  "participant_PCs.csv",
  select=c(1,5:20)
)

pcs[, person_id := trimws(as.character(person_id))]


#related samples
rel <- fread(
  "samples_relatedness_flagged_samples.tsv",
  header=FALSE
)

setnames(rel, "person_id")
rel[, person_id := trimws(as.character(person_id))]


#white British sample list
wb <- fread(
  "aou.v8.classified-rf.wb-classify.report-any.fam",
  header=FALSE
)

wb_ids <- unique(trimws(as.character(wb$V2)))


#physical measures
phys <- fread(
  "physical_measures_hwb.tsv",
  sep="\t"
)

phys[, person_id := trimws(as.character(person_id))]


#labs
labs <- fread(
  "lab_measures_all_v2.tsv",
  sep="\t"
)

labs[, person_id := trimws(as.character(person_id))]


####################################################################################
## Sample QC

related_ids <- unique(rel[!is.na(person_id), person_id])

#remove related individuals from WB sample
keep_ids <- setdiff(wb_ids, related_ids)


#PLINK keep file
plink_keep <- data.table(`#IID`=keep_ids)

fwrite(
  plink_keep,
  "wb_samples_filtered.txt",
  sep="\t"
)


####################################################################################
## Covariates

covs <- merge(
  demo,
  pcs,
  by="person_id",
  all=FALSE
)

#keep unrelated WB samples, male/female
covs <- covs[
  sex %chin% c("Male","Female") &
  person_id %chin% keep_ids
]


covs[, age2 := age^2]

#male=1, female=0
covs[, sex := fifelse(sex == "Male", 1L, 0L)]

covs[, `:=`(
  sex_age=age*sex,
  sex_age2=age2*sex
)]


#center
pc_cols <- names(pcs)[names(pcs)!="person_id"]

cov_names <- c("age",
               "age2",
               "sex",
               "sex_age",
               "sex_age2",
               pc_cols)

covs <- center_covs(covs, cov_names)

covfile <- file.path(
  phdir,
  "aou_covariates.centered.tsv"
)

fwrite(
  covs,
  covfile,
  sep="\t"
)


####################################################################################
## Physical phenotypes

height <- save_pheno(
  phys,
  "height_cm",
  keep_ids,
  file.path(phdir, "aou_height.tsv")
)


weight <- save_pheno(
  phys,
  "weight_kg",
  keep_ids,
  file.path(phdir, "aou_weight.tsv")
)


bmi <- save_pheno(
  phys,
  "bmi",
  keep_ids,
  file.path(phdir, "aou_bmi.tsv")
)


####################################################################################
## Lab phenotypes

#lab file is already one row/person with median values

mono <- save_pheno(
  labs,
  "monocyte_percentage_median",
  keep_ids,
  file.path(phdir, "aou_monocyte_percentage.tsv")
)


baso <- save_pheno(
  labs,
  "basophil_percentage_median",
  keep_ids,
  file.path(phdir, "aou_basophil_percentage.tsv")
)


neut <- save_pheno(
  labs,
  "neutrophil_percentage_median",
  keep_ids,
  file.path(phdir, "aou_neutrophil_percentage.tsv")
)


wbc <- save_pheno(
  labs,
  "white_blood_cell_count_median",
  keep_ids,
  file.path(phdir, "aou_white_blood_cell_count.tsv")
)


rbc <- save_pheno(
  labs,
  "red_blood_cell_count_median",
  keep_ids,
  file.path(phdir, "aou_red_blood_cell_count.tsv")
)


mch <- save_pheno(
  labs,
  "mean_corpuscular_hemoglobin_median",
  keep_ids,
  file.path(phdir, "aou_mean_corpuscular_hemoglobin.tsv")
)


####################################################################################
## Binary phenotypes

asthma <- save_pheno(
  binary,
  "RE_475",
  keep_ids,
  file.path(phdir, "aou_asthma.tsv")
)


t1d <- save_pheno(
  binary,
  "EM_202.1",
  keep_ids,
  file.path(phdir, "aou_t1d.tsv")
)


t2d <- save_pheno(
  binary,
  "EM_202.2",
  keep_ids,
  file.path(phdir, "aou_t2d.tsv")
)


scz <- save_pheno(
  binary,
  "MB_287.1",
  keep_ids,
  file.path(phdir, "aou_schizophrenia.tsv")
)


ad <- save_pheno(
  binary,
  "NS_328.11",
  keep_ids,
  file.path(phdir, "aou_alzheimers.tsv")
)


####################################################################################
## Moments

traits <- list(
  "Height"=height$value,
  "Weight"=weight$value,
  "BMI"=bmi$value,
  "Monocyte percentage"=mono$value,
  "Basophil percentage"=baso$value,
  "Neutrophil percentage"=neut$value,
  "White blood cell count"=wbc$value,
  "Red blood cell count"=rbc$value,
  "Mean corpuscular hemoglobin"=mch$value,
  "Asthma"=asthma$value,
  "Type 1 diabetes"=t1d$value,
  "Type 2 diabetes"=t2d$value,
  "Schizophrenia"=scz$value,
  "Alzheimer's disease"=ad$value
)


moments <- rbindlist(
  lapply(names(traits), function(x) {
    
    m <- calc_moments(traits[[x]])
    
    data.table(
      Trait=x,
      mu0=m["mu0"],
      mu1=m["mu1"],
      mu2=m["mu2"],
      mu3=m["mu3"],
      mu3_std=m["mu3_std"]
    )
    
  }),
  use.names=TRUE
)


fwrite(
  moments,
  file.path(phdir, "aou_moments.tsv"),
  sep="\t"
)


####################################################################################
## Upload

#covariates and sample keep file
system2(
  "gsutil",
  c(
    "-m", "cp",
    covfile,
    "wb_samples_filtered.txt",
    paste0(bucket, "/aou_gwas/")
  )
)


#phenotype files
pheno_files <- c(
  file.path(phdir, "aou_height.tsv"),
  file.path(phdir, "aou_height.plink.tsv"),
  file.path(phdir, "aou_weight.tsv"),
  file.path(phdir, "aou_weight.plink.tsv"),
  file.path(phdir, "aou_bmi.tsv"),
  file.path(phdir, "aou_bmi.plink.tsv"),
  file.path(phdir, "aou_monocyte_percentage.tsv"),
  file.path(phdir, "aou_monocyte_percentage.plink.tsv"),
  file.path(phdir, "aou_basophil_percentage.tsv"),
  file.path(phdir, "aou_basophil_percentage.plink.tsv"),
  file.path(phdir, "aou_neutrophil_percentage.tsv"),
  file.path(phdir, "aou_neutrophil_percentage.plink.tsv"),
  file.path(phdir, "aou_white_blood_cell_count.tsv"),
  file.path(phdir, "aou_white_blood_cell_count.plink.tsv"),
  file.path(phdir, "aou_red_blood_cell_count.tsv"),
  file.path(phdir, "aou_red_blood_cell_count.plink.tsv"),
  file.path(phdir, "aou_mean_corpuscular_hemoglobin.tsv"),
  file.path(phdir, "aou_mean_corpuscular_hemoglobin.plink.tsv"),
  file.path(phdir, "aou_asthma.tsv"),
  file.path(phdir, "aou_asthma.plink.tsv"),
  file.path(phdir, "aou_t1d.tsv"),
  file.path(phdir, "aou_t1d.plink.tsv"),
  file.path(phdir, "aou_t2d.tsv"),
  file.path(phdir, "aou_t2d.plink.tsv"),
  file.path(phdir, "aou_schizophrenia.tsv"),
  file.path(phdir, "aou_schizophrenia.plink.tsv"),
  file.path(phdir, "aou_alzheimers.tsv"),
  file.path(phdir, "aou_alzheimers.plink.tsv"),
  file.path(phdir, "aou_moments.tsv")
)

#bucket
system2(
  "gsutil",
  c(
    "-m", "cp",
    pheno_files,
    paste0(bucket, "/aou_gwas/pheno/")
  )
)