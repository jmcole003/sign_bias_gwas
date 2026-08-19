## Sample QC, split groups for GWAS

library(data.table)
setwd("ukb_group_gwas/phenos/")


####################################################################################
###### Load data

ukb <- fread("UKB.phenos.all.csv")
ukb_qc <- fread("UKB.qc.csv")
ukb <- merge(ukb, ukb_qc, by="eid")

#genotyped samples
gen_samples <- fread("ukb61666_imp_chr1_v3_s487280.sample", sep=" ")
gen_samples <- gen_samples[-1,]

ukb <- ukb[eid %in% gen_samples$ID_1]
setDT(ukb)

####################################################################################
###### Sample QC

#make QC fields numeric
qc_cols <- c("21000-0.0",   # ethnicity
             "21022-0.0",   # age
             "22001-0.0",   # genetic sex
             "22006-0.0",   # genetic ancestry
             "22019-0.0",   # sex chr aneuploidy
             "22020-0.0",   # used in PCA
             "22021-0.0",   # kinship
             "22027-0.0",   # het/missing outlier
             "22028-0.0",   # phasing chr1-22
             "22029-0.0")   # phasing chrX

ukb[, (qc_cols) := lapply(.SD, as.numeric), .SDcols=qc_cols]


#rename PCs
pc_cols<- grep("^22009-0\\.", names(ukb), value=TRUE)

ukb[, (pc_cols) := lapply(.SD, as.numeric), .SDcols=pc_cols]

setnames(ukb,
         old=pc_cols,
         new=paste0("PC",1:length(pc_cols)))


#base QC
base_qc <- (
  !(`22027-0.0` %in% 1) &          #het/missing
  !(`22019-0.0` %in% 1) &          #sex chr aneuploidy
  !(`22021-0.0` %in% c(-1,10)) &   #kinship
  (`22020-0.0` %in% 1)             #used in PCA
)


####################################################################################
###### Ancestry

#white British reference
wb_idx<- which(
  ukb$`22006-0.0` == 1 &
  base_qc
)

pc_ellipse <- paste0("PC",1:6)

pc_means <- sapply(pc_ellipse, function(x)
  mean(ukb[wb_idx, get(x)], na.rm=TRUE))

pc_sds <- sapply(pc_ellipse, function(x)
  sd(ukb[wb_idx, get(x)], na.rm=TRUE))


#distance from WB center
dd <- rep(0,nrow(ukb))

for(i in 1:6){
  pc <- ukb[[pc_ellipse[i]]]
  dd <- dd + (pc-pc_means[i])^2/(pc_sds[i]^2)
}

ukb[, dd_brit := dd]

#self reported white/missing ethnicity
ethnic <- ukb$`21000-0.0`

white <- ethnic %in% c(1,1001,1002,1003)
ethnic_miss <- ethnic %in% c(-3,-1) | is.na(ethnic)

ancestry_keep <- ukb$dd_brit < 7^2 & (white | ethnic_miss)


####################################################################################
###### Final GWAS sample

ukb_final <- ukb[
  base_qc &
  (`22028-0.0` %in% 1) &
  (`22029-0.0` %in% 1) &
  !is.na(`21022-0.0`) &
  ancestry_keep
]


#sample keep file
keep_qc <- ukb_final[, .(
  FID=eid,
  IID=eid
)]

fwrite(keep_qc,
       "ukb_qc_samples.txt",
       sep="\t")


####################################################################################
###### Covariates

ukb_final[, age := as.numeric(`21022-0.0`)]
ukb_final[, sex := as.numeric(`22001-0.0`)]

ukb_final <- ukb_final[
  !is.na(age) &
  !is.na(sex)
]


#age terms
ukb_final[, age2 := age^2]
ukb_final[, sex_age := sex*age]
ukb_final[, sex_age2 := sex*age2]


#covariate file
pc_use <- paste0("PC",1:20)

covars <- ukb_final[, c(
  list(FID=eid,
       IID=eid,
       sex=sex,
       age=age,
       age2=age2,
       sex_age=sex_age,
       sex_age2=sex_age2),
  mget(pc_use)
)]

fwrite(covars,
       "ukb.covariates.txt",
       sep="\t")


####################################################################################
###### Binary phenotypes

#main ICD10 fields
icd_cols <- grep("^41202-", names(ukb_final), value=TRUE)

icd <- melt(
  ukb_final,
  id.vars="eid",
  measure.vars=icd_cols,
  value.name="icd10",
  variable.factor=FALSE
)

icd <- icd[!is.na(icd10) & icd10!=""]

icd[, icd10 := toupper(trimws(as.character(icd10)))]
icd[, icd3 := substr(icd10,1,3)]


#case IDs
ad_ids <- icd[icd3=="G30", unique(eid)]
ast_ids <- icd[icd3=="J45", unique(eid)]

t1d_ids <- icd[icd3=="E10", unique(eid)]
t2d_ids <- icd[icd3=="E11", unique(eid)]

#T1D gets priority
t2d_ids <- setdiff(t2d_ids,t1d_ids)

scz_ids <- icd[icd3=="F20", unique(eid)]


#binary phenos
ukb_final[, AD := as.integer(eid %in% ad_ids)]
ukb_final[, ASTHMA := as.integer(eid %in% ast_ids)]
ukb_final[, T1D := as.integer(eid %in% t1d_ids)]
ukb_final[, T2D := as.integer(eid %in% t2d_ids)]
ukb_final[, SCZ := as.integer(eid %in% scz_ids)]


#case counts
summary_counts <- data.table(
  trait=c("AD","T1D","T2D","SCZ","ASTHMA"),
  
  n_cases=c(
    sum(ukb_final$AD),
    sum(ukb_final$T1D),
    sum(ukb_final$T2D),
    sum(ukb_final$SCZ),
    sum(ukb_final$ASTHMA)
  )
)

summary_counts[, n_total := nrow(ukb_final)]
summary_counts[, n_controls := n_total-n_cases]

summary_counts[, neale_cases :=
                 c(119,583,888,198,1693)]

summary_counts[, neale_controls :=
                 c(361075,360611,360306,360996,359501)]

print(summary_counts)


####################################################################################
###### Quantitative phenotypes

ukb_final[, height := as.numeric(`50-0.0`)]
ukb_final[, bmi := as.numeric(`21001-0.0`)]
ukb_final[, weight := as.numeric(`21002-0.0`)]

ukb_final[, mono_pct := as.numeric(`30190-0.0`)]
ukb_final[, neutro_pct := as.numeric(`30200-0.0`)]
ukb_final[, baso_pct := as.numeric(`30220-0.0`)]

ukb_final[, wbc := as.numeric(`30000-0.0`)]
ukb_final[, rbc := as.numeric(`30010-0.0`)]
ukb_final[, mch := as.numeric(`30050-0.0`)]


####################################################################################
###### Phenotype file

#plink binary coding: 2=control, 3=case
ukb_final[, AD := AD+2L]
ukb_final[, T1D := T1D+2L]
ukb_final[, T2D := T2D+2L]
ukb_final[, SCZ := SCZ+2L]
ukb_final[, ASTHMA := ASTHMA+2L]


pheno_plink <- ukb_final[, .(
  FID=eid,
  IID=eid,
  
  AD,
  T1D,
  T2D,
  SCZ,
  ASTHMA,
  
  height,
  bmi,
  weight,
  
  mono_pct,
  neutro_pct,
  baso_pct,
  
  wbc,
  rbc,
  mch
)]


fwrite(pheno_plink,
       "ukb_phenotypes.txt",
       sep="\t",
       na="NA")


####################################################################################
###### Split QC samples in half

qc <- copy(keep_qc)

qc <- qc[sample(.N)]

n <- nrow(qc)
n1 <- floor(n/2)

qc_group1 <- qc[1:n1]
qc_group2 <- qc[(n1+1):n]


fwrite(qc_group1,
       "ukb_samples_group1.txt",
       sep="\t")

fwrite(qc_group2,
       "ukb_samples_group2.txt",
       sep="\t")


####################################################################################
###### Split QC samples into thirds

qc <- copy(keep_qc)

qc <- qc[sample(.N)]

qc[, group := rep(1:3,length.out=.N)]


qc_small1 <- qc[group==1, .(FID,IID)]
qc_small2 <- qc[group==2, .(FID,IID)]
qc_small3 <- qc[group==3, .(FID,IID)]


fwrite(qc_small1,
       "ukb_samples_smallgroup1.txt",
       sep="\t")

fwrite(qc_small2,
       "ukb_samples_smallgroup2.txt",
       sep="\t")

fwrite(qc_small3,
       "ukb_samples_smallgroup3.txt",
       sep="\t")


####################################################################################
###### Check binary counts in splits

binary_traits <- c("AD","T1D","T2D","SCZ","ASTHMA")


split_counts <- function(ids,name){
  
  x <- ukb_final[eid %in% ids]
  
  out <- rbindlist(lapply(binary_traits,function(tr){
    
    v <- x[[tr]]
    
    data.table(
      split=name,
      trait=tr,
      n_total=length(v),
      n_cases=sum(v==3,na.rm=TRUE),
      n_controls=sum(v==2,na.rm=TRUE),
      n_missing=sum(is.na(v))
    )
    
  }))
  
  return(out)
}


summary_splits <- rbindlist(list(
  
  split_counts(qc_group1$IID,"big_group1"),
  split_counts(qc_group2$IID,"big_group2"),
  
  split_counts(qc_small1$IID,"small_group1"),
  split_counts(qc_small2$IID,"small_group2"),
  split_counts(qc_small3$IID,"small_group3")
  
))


print(summary_splits)