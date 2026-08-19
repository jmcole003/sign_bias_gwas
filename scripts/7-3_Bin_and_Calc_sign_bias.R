#!/usr/bin/env Rscript

######## Binning/threshold
library(data.table)

#arguments
args <- commandArgs(trailingOnly=TRUE)

file_arg <- args[which(args=="--file")+1]
type_arg <- args[which(args=="--type")+1]
method_arg <- args[which(args=="--method")+1]
mode_arg <- args[which(args=="--mode")+1]

low_idx <- which(args=="--low")
high_idx <- which(args=="--high")

T_low <- if(length(low_idx)>0) as.numeric(args[low_idx+1]) else NA
T_high <- if(length(high_idx)>0) as.numeric(args[high_idx+1]) else NA

trait_name <- sub(".*\\.(.*?)\\..*","\\1",basename(file_arg))


####################################################################################
#bins

bin_breaks_maf <- c(0,5e-4,1e-3,0.002,0.003,0.005,0.01,
                    0.02,0.05,0.1,0.3,0.5)

bin_breaks_daf <- c(0,5e-4,1e-3,0.002,0.003,0.005,0.01,
                    0.02,0.05,0.1,0.3,0.5,0.7,0.9,
                    0.95,0.98,0.99,0.995,0.998,1)


#read
if(type_arg=="MAF"){
  dat <- fread(file_arg,select=c("block_id","MAF","eta","pval"))
  dat$af <- dat$MAF
}

if(type_arg=="DAF"){
  dat <- fread(file_arg,select=c("block_id","DAF","eta","pval"))
  dat$af <- dat$DAF
}


####################################################################################
###### bins

if(mode_arg=="bin"){

  if(type_arg=="MAF") breaks <- bin_breaks_maf
  if(type_arg=="DAF") breaks <- bin_breaks_daf

  bin_labels <- paste0("bin",1:(length(breaks)-1))

  dat$bin_id <- cut(dat$af,
                    breaks=breaks,
                    labels=bin_labels,
                    include.lowest=TRUE)

  res <- vector("list",length(bin_labels))


  for(i in 1:length(bin_labels)){

    x <- dat[dat$bin_id==bin_labels[i],]

    #blocks with >=10 SNPs
    cnt <- x[,.N,by=block_id]
    blocks <- cnt[N>=10]$block_id
    x <- x[x$block_id %in% blocks,]


    if(nrow(x)==0){

      res[[i]] <- data.frame(
        bin=bin_labels[i],
        bin_range=paste0("[",breaks[i],",",breaks[i+1],")"),
        n_snps=0,
        n_blocks=0,
        mean_eta=NA_real_,
        se=NA_real_,
        lower_ci=NA_real_,
        upper_ci=NA_real_
      )

      next
    }


    n_snps <- nrow(x)
    blocks <- unique(x$block_id)
    n_blocks <- length(blocks)


    #random
    if(method_arg=="random"){

      block_eta<- split(x$eta,x$block_id)
      eta_rep <- numeric(1000)

      for(r in 1:1000){

        eta_use <- sapply(block_eta,function(z) sample(z,1))

        eta_rep[r] <- sum(eta_use)/sum(abs(eta_use))
      }

      mean_eta <- mean(eta_rep)
      se <- sd(eta_rep)
    }


    #lowest p
    if(method_arg=="sig"){

      block_dat <- split(x,x$block_id)

      eta_use <- sapply(block_dat,function(z){
        z$eta[which.min(z$pval)]
      })

      mean_eta <- sum(eta_use)/sum(abs(eta_use))


      #bootstrap blocks
      eta_rep <- numeric(1000)

      for(r in 1:1000){

        bs <- sample(blocks,length(blocks),replace=TRUE)
        eta_use <- c()

        for(b in bs){

          z<- block_dat[[as.character(b)]]

          eta_use <- c(
            eta_use,
            z$eta[which.min(z$pval)]
          )
        }

        eta_rep[r] <- sum(eta_use)/sum(abs(eta_use))
      }

      se <- sd(eta_rep,na.rm=TRUE)
    }


    res[[i]] <- data.frame(
      bin=bin_labels[i],
      bin_range=paste0("[",breaks[i],",",breaks[i+1],")"),
      n_snps=n_snps,
      n_blocks=n_blocks,
      mean_eta=mean_eta,
      se=se,
      lower_ci=mean_eta-1.96*se,
      upper_ci=mean_eta+1.96*se
    )
  }


  bin_df<- do.call(rbind,res)
  rownames(bin_df) <- NULL

  bin_df$trait <- trait_name
  bin_df$type <- type_arg
  bin_df$method <- method_arg


  out_file <- paste0("bin_",trait_name,"_",type_arg,"_",
                     method_arg,".txt")

  fwrite(bin_df,out_file,sep="\t",quote=FALSE,row.names=FALSE)

}


####################################################################################
###### threshold

if(mode_arg=="threshold"){

  low <- if(is.na(T_low)) 0 else T_low


  #filter
  if(type_arg=="MAF"){
    x <- dat[dat$af<low,]
  }

  if(type_arg=="DAF"){

    if(!is.na(T_high)){
      x <- dat[dat$af<low | dat$af>T_high,]
    } else {
      x <- dat[dat$af<low,]
    }
  }


  #blocks with >=10 SNPs
  cnt <- x[,.N,by=block_id]
  blocks <- cnt[N>=10]$block_id
  x <- x[x$block_id %in% blocks,]


  if(nrow(x)==0){

    out_df <- data.frame(
      mean_eta=NA_real_,
      se=NA_real_,
      lower_ci=NA_real_,
      upper_ci=NA_real_,
      n_snps=0,
      n_blocks=0
    )

  } else {

    n_snps <- nrow(x)
    blocks <- unique(x$block_id)
    n_blocks <- length(blocks)


    #random
    if(method_arg=="random"){

      block_eta <- split(x$eta,x$block_id)
      eta_rep <- numeric(1000)

      for(r in 1:1000){

        eta_use <- sapply(block_eta,function(z) sample(z,1))

        eta_rep[r] <- sum(eta_use)/sum(abs(eta_use))
      }

      mean_eta <- mean(eta_rep)
      se <- sd(eta_rep)
    }


    #lowest p
    if(method_arg=="sig"){

      block_dat <- split(x,x$block_id)

      eta_use <- sapply(block_dat,function(z){
        z$eta[which.min(z$pval)]
      })

      mean_eta <- sum(eta_use)/sum(abs(eta_use))


      #bootstrap
      eta_rep <- numeric(1000)

      for(r in 1:1000){

        bs <- sample(blocks,length(blocks),replace=TRUE)
        eta_use <- c()

        for(b in bs){

          z <- block_dat[[as.character(b)]]

          eta_use <- c(
            eta_use,
            z$eta[which.min(z$pval)]
          )
        }

        eta_rep[r] <- sum(eta_use)/sum(abs(eta_use))
      }

      se <- sd(eta_rep,na.rm=TRUE)
    }


    out_df <- data.frame(
      mean_eta=mean_eta,
      se=se,
      lower_ci=mean_eta-1.96*se,
      upper_ci=mean_eta+1.96*se,
      n_snps=n_snps,
      n_blocks=n_blocks
    )
  }


  out_df$trait <- trait_name
  out_df$type <- type_arg
  out_df$method <- method_arg
  out_df$T_low <- T_low
  out_df$T_high <- T_high


  low_str <- if(!is.na(T_low)) paste0("_",T_low) else ""
  high_str <- if(!is.na(T_high)) paste0("_",T_high) else ""

  out_file <- paste0("threshold_",trait_name,"_",type_arg,"_",
                     method_arg,low_str,high_str,".txt")


  fwrite(out_df,out_file,sep="\t",quote=FALSE, row.names=FALSE)

}