## Simulation of skewed cohorts (scheme A)

# n_pop: population size
# n_pairs: number of matched pairs (increasing/decreasing)
# h2_target: target SNP h2
# n_samp: cohort size
# tau: mixture weight for bin sampling
# gamma: concentration near q
# block_size: pseudo-ld block size
# maf_max: MAF threshold
# skew_tol: tolerance level for target skew
# max_tries: maximum sampling attempts for each rep

library(dplyr)
library(tibble)
library(ggplot2)
library(ashr)
library(cowplot)

# params
n_pop<-2000000 
n_pairs<-2000 
h2_target<-0.6 
n_bins<-200
targets<-1:10
n_rep<-20
n_samp<-10000 
tau<-0.01 
gamma<-20 
q_grid<-seq(0, 0.95, by=0.0001)
#q_grid<-seq(0, 0.99, by=0.0001)
#q_grid<-seq(0, 0.90, by=0.0001)
tune_seed<-23276
base_seed<-1000
gwas_chunk<-200
block_size<-6 
ci_level<-0.95
maf_max<-0.01 
skew_tol<-0.5 
max_tries<-50 

skew3<-function(x){
  m<-mean(x)
  s<-sd(x)
  if (s==0){
    sk<-0
  } else {
    sk<-mean(((x-m)/s)^3)
  }
  return(sk)
}

### Build phenotype and genotypes
build_population<-function(N_pop=140000,
                           M_pairs=1000,
                           maf_min=5e-4,
                           maf_max=8e-3,
                           h2=0.6,
                           chunk=200,
                           seed=1) {
  set.seed(seed)
  
  # number of snps and mafs
  M<-2*M_pairs
  p_pair<-10^(runif(M_pairs,
                    log10(maf_min), 
                    log10(maf_max)))
  p<-rep(p_pair, each=2)
  
  sum_varG<-sum(2*p*(1-p))
  bmag<-sqrt(h2/sum_varG)
  
  beta <- rep(c(+bmag,-bmag), times=M_pairs) #equal increasing/decreasing
  
  G_raw<-matrix(as.raw(0), nrow=N_pop, ncol=M) # genotypes
  #dim(G_raw) 
  gval<-numeric(N_pop)
  
  for(startsnp in seq(1,M,by=chunk)){
    snps <- startsnp:min(startsnp+chunk-1, M)
    G <- matrix(rbinom(N_pop*length(snps),
                       2,
                       rep(p[snps], each=N_pop)),
              nrow=N_pop,ncol=length(snps))
    
    G_raw[,snps]<-as.raw(G)
    gval <- gval + as.vector(G%*%beta[snps])
  }
  
  # phenotype
  # dim(gval)
  e<-rnorm(N_pop, sd=sqrt(1-h2))
  Y_pop <- gval + e
  Y_pop <- as.numeric(Y_pop)
  
  list(
    N_pop=N_pop,M=M,p=p,beta=beta,
    G_raw=G_raw,Y_pop=Y_pop,h2=h2,bmag = bmag
  )
}

##make phenotype bins
prep_bins <- function(Y,B=200){
  
  probs <- seq(0,1,length.out=B+1)
  breaks <- quantile(Y,probs)
  breaks <- unique(as.numeric(breaks))
  
  n_bins <- length(breaks)-1
  
  #binning 
  bin <- cut(Y,breaks=breaks,include.lowest=TRUE,labels=FALSE)
  inds <- split(1:length(Y),bin)
  sizes <- sapply(inds, length)
  
  #midpoint
  r_mid <- ((1:n_bins)-0.5)/n_bins
  
  out<-list("Y"=Y,
            "brks"=breaks,
            "idx_by_bin"=inds,
            "sizes"=sizes,
            "r_mid"=r_mid,
            "B"=n_bins)
  
  return(out)
}

# allocate number drawn from bins
allocate_counts <- function(N_samp, prob, sizes){
  
  prob <- prob/sum(prob)
  counts <- as.vector(rmultinom(1,N_samp,prob=prob))
  over <- which(counts>sizes)

  while(length(over)>0){
    extra <- sum(counts[over]-sizes[over])
    counts[over] <- sizes[over]
    space <- sizes-counts

    #redistribute extra individuals
    prob2 <- space*prob
    prob2 <- prob2/sum(prob2)

    add <- as.vector(rmultinom(1,extra,prob=prob2))
    counts <- counts+add

    over <- which(counts>sizes)
  }

  return(counts)
}


### Cohort sampler 
sample_cohort<-function(prep,
                        N_samp=7000,
                        q=0,
                        gamma=5,
                        tau=0.05,
                        seed=1) {
  set.seed(seed)
  
  r<-prep$r_mid
  keep_bins<-which(r >= q)
  r_keep<-(r[keep_bins]-q)/(1-q+1e-12)
  
  # weights high near cutoff, low near top
  w<-(1-r_keep)^gamma
  p_keep<-w/sum(w)
  
  K<-length(keep_bins)
  prob<-rep(0, prep$B)
  prob[keep_bins]<-tau*rep(1/K, K)+(1-tau)*p_keep
  
  cnt<-allocate_counts(N_samp, prob, prep$sizes)
  
  idx<-unlist(mapply(function(v, k) if (k > 0) sample(v, k, FALSE) else integer(0),
                     prep$idx_by_bin, cnt,
                     SIMPLIFY=FALSE, USE.NAMES=FALSE))
  y_samp<-prep$Y[idx]
  
  list(idx=idx, Y_samp=y_samp,
       mean=mean(y_samp), skew=skew3(y_samp),
       q=q, gamma=gamma,tau=tau)
}

# GWAS (linreg)
gwas_all_snps <- function(pop,idx,chunk=200){
  
  n <- length(idx)
  y <- pop$Y_pop[idx]
  y <- y-mean(y)
  Syy <- sum(y^2)
  
  n_snps <- ncol(pop$G_raw)
  
  #preallocate vectors
  betahat <- rep(NA,n_snps)
  se <- rep(NA,n_snps)
  af <- rep(NA,n_snps)
  maf <- rep(NA,n_snps)
  mac <- rep(NA,n_snps)
  
  #run SNPs in chunks
  for(start_snp in seq(1,n_snps,by=chunk)){
    
    end_snp <- min(start_snp+chunk-1,n_snps)
    snps <- start_snp:end_snp
  
    G <- matrix(as.integer(pop$G_raw[idx,snps,drop=FALSE]),
                nrow=n,ncol=length(snps))
    
    allele_count <- colSums(G)
    
    af_snp <- allele_count/(2*n)
    maf_snp <- pmin(af_snp,1-af_snp)
    mac_snp <- pmin(allele_count,2*n-allele_count)
    
    #allele freqs
    af[snps] <- af_snp
    maf[snps] <- maf_snp
    mac[snps] <- mac_snp
    
    keep <- mac_snp>=1 & maf_snp>0
    if(!any(keep)) next
    
    #linear
    gmean <- allele_count/n
    Sxx <- colSums(G^2)-n*(gmean^2)
    Sxy <- as.numeric(crossprod(G,y))
    
    b <- Sxy/Sxx
    
    #orient effect for minor allele
    flip <- af_snp>0.5
    b[flip] <- -b[flip]
    

    SSE <- Syy-(Sxy^2)/Sxx
    sigma2 <- SSE/(n-2)
    s <- sqrt(sigma2/Sxx)
    
    betahat[snps] <- b
    se[snps] <- s
  }
  
  #z score
  z <- betahat/se
  
  #out
  out<-data.frame("snp"=1:n_snps,
                  "af"=af,
                  "maf"=maf,
                  "mac"=mac,
                  "betahat"=betahat,
                  "se"=se,
                  "z"=z)
  
  return(out)
}


# ASH derived sign bias
run_ash_sb<-function(bhat,se){
  sb<-rep(NA,length(bhat))
  ok<-is.finite(bhat)&is.finite(se)&(se>0)
  
  #run ash
  fit<-ash(betahat=bhat[ok],sebetahat=se[ok],method="fdr",mixcompdist="normal",optmethod="mixSQP")
  
  if(!is.null(fit$result)&&all(c("PositiveProb","NegativeProb") %in% colnames(fit$result))){
    sb_ok<-as.numeric(fit$result$PositiveProb-fit$result$NegativeProb)
  } else {
    p_ge0<-as.numeric(get_posterior_prob(fit,l=0,u=Inf))
    p_le0<-as.numeric(get_posterior_prob(fit,l=-Inf,u=0))
    sb_ok<-as.numeric(p_ge0-p_le0)
  }
  
  sb[ok]<-sb_ok
  return(sb)
}

#cohort sign bias
# S=sum(sign score)/sum(abs(sign score))
cohort_sign_bias_from_sb <- function(sb){
  use <- is.finite(sb)
  sb <- sb[use]
  
  #calc total
  if(length(sb)==0) return(NA)
  total <- sum(abs(sb),na.rm=TRUE)
  
  #run check 2 / calc bias
  if(total==0) return(NA)
  bias <- sum(sb, na.rm=TRUE)/total
  
  return(bias)
}

#true effect sign, oriented to minor allele
true_minor_sign <- function(allelefreq, beta){
  true_sign <- sign(beta)
  true_sign[allelefreq>0.5] = -true_sign[allelefreq>0.5]
  return(true_sign)
}

#lowest p SNP from each pseudo ld block
min_p_per_block <- function(gwas, inds, 
                            block_size=6){
  
  #p <- 2*pnorm(testGWAS$Z)
  p <- 2*pnorm(-abs(gwas$z[inds])) #gwas p
  
  block <- ((gwas$snp[inds]-1)%/%block_size)+1
  blocks <- sort(unique(block))
  
  winners <- rep(NA,length(blocks))
  
  #lowest p SNP 
  for(i in 1:length(blocks)){
    in_block <- which(block==blocks[i])
    min_p <- which.min(p[in_block])
    winners[i] <- inds[in_block[min_p]]
  }
  
  return(winners)
}

#sample until skew is within tolerance
sample_to_skew <- function(prep,N_samp,q,gamma,tau,seed,
                           target_sk,skew_tol,max_tries,
                           give_sample=NULL){
  
  best_err<-Inf
  
  if(!is.null(give_sample)){
    
    sk<-skew3(give_sample$Y_samp)
    err<-abs(sk-target_sk)
    
    best<-give_sample
    best_err<-err
    
    if(err<=skew_tol){
      return(list(samp=give_sample,skew=sk,err=err,
                  in_tol=TRUE,n_tries=1))
    }
    
    first_try<-2
    
  } else {
    
    first_try<-1
  }
  
  #resample
  for(i in first_try:max_tries){
    
    samp<-sample_cohort(prep,N_samp=N_samp,q=q,gamma=gamma,tau=tau,
                        seed=seed+1000000*(i-1))
    
    sk<-skew3(samp$Y_samp)
    err<-abs(sk-target_sk)
    
    if(err<best_err){
      best<-samp
      best_err<-err
    }
    
    if(err<=skew_tol){
      return(list(samp=samp,skew=sk,err=err,
                  in_tol=TRUE,n_tries=i))
    }
  }
  
  out<-list(samp=best,
            skew=skew3(best$Y_samp),
            err=best_err,
            in_tol=FALSE,
            n_tries=max_tries)
  
  return(out)
}


#fit ASH and calculate sign bias for one cohort
one_cohort_bias <- function(pop,prep,target_sk,q,gamma,N_samp,
                            tau=0.05,seed=1,gwas_chunk=200,
                            block_size=6,maf_max=Inf,
                            skew_tol=Inf,max_tries=1){
  
  #sample cohort
  pick <- sample_to_skew(prep,N_samp=N_samp,q=q,gamma=gamma,tau=tau,
                         seed=seed,target_sk=target_sk,
                         skew_tol=skew_tol,max_tries=max_tries)
  
  samp <- pick$samp
  
  #if target skew was not reached
  if(!pick$in_tol){
    
    out<-data.frame(target_sk=target_sk,
                    skew_obs=pick$skew,
                    skew_err=pick$err,
                    in_tol=FALSE,
                    n_tries=pick$n_tries,
                    q=q,
                    gamma=gamma,
                    bias_sig=NA,
                    true_bias_sample=NA,
                    n_used=0,
                    seed=seed)
    
    return(out)
  }
  
  #GWAS
  gwas <- gwas_all_snps(pop,idx=samp$idx,chunk=gwas_chunk)
  
  #SNPs that can be used for ASH
  use_ash <- !is.na(gwas$betahat) & !is.na(gwas$se) &
    gwas$se>0 & gwas$mac>=1
  
  sb <- rep(NA,nrow(gwas))
  sb[use_ash] <- run_ash_sb(gwas$betahat[use_ash],gwas$se[use_ash])
  
  #rare SNPs used for sign bias
  use <- use_ash & gwas$maf>0 & gwas$maf<=maf_max &
    !is.na(sb) & !is.na(gwas$z)
  
  inds <- which(use)
  
  #true sign bias
  beta_true <- pop$beta[gwas$snp[inds]]
  true_sign <- true_minor_sign(gwas$af[inds],beta_true)
  true_bias <- mean(true_sign)
  
  #lowest p SNP from each pseudo-block
  winners <- min_p_per_block(gwas,inds,block_size=block_size)
  
  #estimated sign bias
  bias_sig <- cohort_sign_bias_from_sb(sb[winners])
  
  #out
  out<-data.frame(target_sk=target_sk,
                  skew_obs=pick$skew,
                  skew_err=pick$err,
                  in_tol=TRUE,
                  n_tries=pick$n_tries,
                  q=q,
                  gamma=gamma,
                  bias_sig=bias_sig,
                  true_bias_sample=true_bias,
                  n_used=length(winners),
                  seed=seed)
  
  return(out)
}

##find q that gives closest target skew
tune_q_for_skew <- function(prep,target_sk,N_samp,gamma,tau,
                            q_grid=seq(0,0.95,by=0.01),
                            seed=1,skew_tol=Inf){
  
  best_err <- Inf
  best_q <- NA
  best_skew <- NA
  
  for(q in q_grid){
    
    samp <- sample_cohort(prep,N_samp=N_samp,q=q,
                          gamma=gamma,tau=tau,seed=seed)
    
    sk <- skew3(samp$Y_samp)
    err <- abs(sk-target_sk)
    
    if(err<best_err){
      best_err <- err
      best_q <- q
      best_skew <- sk
    }
  }
  
  in_tol <- best_err<=skew_tol
  
  if(!in_tol){
    best_q <- NA
  }
  
  #out
  out<-data.frame(target_sk=target_sk,
                  q=best_q,
                  gamma=gamma,
                  tau=tau,
                  skew_pilot=best_skew,
                  err=best_err,
                  in_tol=in_tol)
  
  return(out)
}

##run simulations for each target skew
simulate_signbias_vs_skew <- function(pop,prep,targets=1:10,n_rep=20,
                                      N_samp=10000,tau=0.01,gamma=20,
                                      tune_seed=1,q_grid=seq(0,0.95,by=0.01),
                                      gwas_chunk=200,block_size=6,
                                      maf_max=Inf,ci_level=0.95,
                                      base_seed=1000,skew_tol=Inf,
                                      max_tries=1){
  
  #### tune q
  
  tuning_list <- list()
  
  for(i in 1:length(targets)){
    
    tuning_list[[i]] <- tune_q_for_skew(
      prep,
      target_sk=targets[i],
      N_samp=N_samp,
      gamma=gamma,
      tau=tau,
      q_grid=q_grid,
      seed=tune_seed+i,
      skew_tol=skew_tol
    )
  }
  
  tuning <- do.call(rbind,tuning_list)
  
  
  #### run replicates
  
  reps_list <- list()
  k <- 1
  
  for(i in 1:length(targets)){
    
    target <- targets[i]
    q <- tuning$q[tuning$target_sk==target]
    
    #skip target if tuning did not get within tolerance
    if(is.na(q)) next
    
    for(r in 1:n_rep){
      
      reps_list[[k]] <- one_cohort_bias(
        pop,prep,
        target_sk=target,
        q=q,
        gamma=gamma,
        N_samp=N_samp,
        tau=tau,
        seed=base_seed+10000*i+r,
        gwas_chunk=gwas_chunk,
        block_size=block_size,
        maf_max=maf_max,
        skew_tol=skew_tol,
        max_tries=max_tries
      )
      
      k <- k+1
    }
  }
  
  reps <- do.call(rbind,reps_list)
  
  
  #### summarize
  
  alpha <- 1-ci_level
  sum_list <- list()
  
  for(i in 1:length(targets)){
    
    target <- targets[i]
    x <- reps[reps$target_sk==target,]
    
    #nothing was run for this target
    if(nrow(x)==0) next
    
    n_ok <- sum(!is.na(x$bias_sig))
    
    mean_bias <- mean(x$bias_sig,na.rm=TRUE)
    sd_bias <- sd(x$bias_sig,na.rm=TRUE)
    se_bias <- sd_bias/sqrt(n_ok)
    
    mean_true <- mean(x$true_bias_sample,na.rm=TRUE)
    sd_true <- sd(x$true_bias_sample,na.rm=TRUE)
    se_true <- sd_true/sqrt(n_ok)
    
    tcrit <- qt(1-alpha/2,df=n_ok-1)
    
    sum_list[[i]] <- data.frame(
      target_sk=target,
      n_rep=n_ok,
      n_in_tol=sum(x$in_tol),
      mean_skew_obs=mean(x$skew_obs,na.rm=TRUE),
      
      mean_bias=mean_bias,
      sd_bias=sd_bias,
      se_bias=se_bias,
      ci_low_bias=mean_bias-tcrit*se_bias,
      ci_high_bias=mean_bias+tcrit*se_bias,
      
      mean_true=mean_true,
      sd_true=sd_true,
      se_true=se_true,
      ci_low_true=mean_true-tcrit*se_true,
      ci_high_true=mean_true+tcrit*se_true
    )
  }
  
  summ <- do.call(rbind,sum_list)
  
  #add q used for each target
  q_used <- tuning$q[match(summ$target_sk,tuning$target_sk)]
  skew_pilot <- tuning$skew_pilot[match(summ$target_sk,tuning$target_sk)]
  
  summ$q <- q_used
  summ$gamma <- gamma
  summ$tau <- tau
  summ$skew_pilot <- skew_pilot
  
  
  #out
  out<-list("tuning"=tuning,
            "reps"=reps,
            "summary"=summ)
  
  return(out)
}


####################################################################################
## RUN

#simulate population
pop<-build_population(N_pop=n_pop,M_pairs=n_pairs,
                      maf_max=0.1,h2=h2_target,seed=1)

#make quantile bins
prep<-prep_bins(pop$Y_pop,B=n_bins)

#run simulations
out<-simulate_signbias_vs_skew(
  pop,prep,
  targets=targets,
  n_rep=n_rep,
  N_samp=n_samp,
  tau=tau,
  gamma=gamma,
  q_grid=q_grid,
  tune_seed=tune_seed,
  gwas_chunk=gwas_chunk,
  block_size=block_size,
  maf_max=maf_max,
  ci_level=ci_level,
  base_seed=base_seed,
  skew_tol=skew_tol,
  max_tries=max_tries
)

####################################################################################
#add baseline skew=0

reps0<-lapply(1:n_rep,function(r){
  
  one_cohort_bias(
    pop,prep,
    target_sk=0,
    q=0,
    gamma=0,
    N_samp=n_samp,
    tau=1,
    seed=900000+r,
    gwas_chunk=gwas_chunk,
    block_size=block_size,
    maf_max=maf_max,
    skew_tol=skew_tol,
    max_tries=max_tries
  )
})

reps0<-do.call(rbind,reps0)

#combine baseline with skewed cohorts
reps_all<-rbind(reps0,out$reps)

alpha<-1-ci_level

#example skewed cohort for plotting
example_cohort<-sample_cohort(prep,N_samp=n_samp,q=0.46,
                              gamma=gamma,tau=tau,seed=199)

####################################################################################
##### Plotting


## Population distribution

dens_pop<-density(pop$Y_pop)
pop_plot<-data.frame(x=dens_pop$x,y=dens_pop$y)

pA<-ggplot(pop_plot,aes(x,y))+
  geom_area(fill="#FDE992",color=NA)+
  geom_line(color="black",linewidth=0.6)+
  geom_vline(xintercept=-0.1,linewidth=0.8,linetype=2)+
  theme_classic()+
  theme(axis.title=element_blank(),
        axis.text=element_text(size=10))


## Example cohort distribution

y_example<-example_cohort$Y_samp
dens_cohort<-density(y_example,adjust=3)
cohort_plot<-data.frame(x=dens_cohort$x,y=dens_cohort$y)

pB<-ggplot(cohort_plot,aes(x,y))+
  geom_area(fill="#EFEFEF",color=NA)+
  geom_line(color="black",linewidth=0.6)+
  theme_classic()+
  theme(axis.title=element_blank(),
        axis.text=element_text(size=10))


####################################################################################
## Sign bias

res_sim<-out$summary


## skew=0 baseline

x0<-reps_all[reps_all$target_sk==0,]

n0<-sum(!is.na(x0$bias_sig))

mean_bias<-mean(x0$bias_sig,na.rm=TRUE)
sd_bias<-sd(x0$bias_sig,na.rm=TRUE)
se_bias<-sd_bias/sqrt(n0)

mean_true<-mean(x0$true_bias_sample,na.rm=TRUE)
sd_true<-sd(x0$true_bias_sample,na.rm=TRUE)
se_true<-sd_true/sqrt(n0)

tcrit<-qt(1-alpha/2,df=n0-1)

skew0<-data.frame(
  target_sk=0,
  n_rep=n0,
  n_in_tol=sum(x0$in_tol),
  mean_skew_obs=mean(x0$skew_obs,na.rm=TRUE),
  
  mean_bias=mean_bias,
  sd_bias=sd_bias,
  se_bias=se_bias,
  ci_low_bias=mean_bias-tcrit*se_bias,
  ci_high_bias=mean_bias+tcrit*se_bias,
  
  mean_true=mean_true,
  sd_true=sd_true,
  se_true=se_true,
  ci_low_true=mean_true-tcrit*se_true,
  ci_high_true=mean_true+tcrit*se_true
)


## add baseline

#only need columns used below
plot_res<-rbind(
  res_sim[,c("target_sk","mean_skew_obs","mean_bias",
             "ci_low_bias","ci_high_bias","mean_true",
             "ci_low_true","ci_high_true")],
  
  skew0[,c("target_sk","mean_skew_obs","mean_bias",
           "ci_low_bias","ci_high_bias","mean_true",
           "ci_low_true","ci_high_true")]
)

plot_res<-plot_res[order(plot_res$target_sk),]


## data for each line

pop_line<-data.frame(
  x=plot_res$mean_skew_obs,
  mean=0
)

true_line<-data.frame(
  x=plot_res$mean_skew_obs,
  mean=plot_res$mean_true,
  lower=plot_res$ci_low_true,
  upper=plot_res$ci_high_true
)

est_line<-data.frame(
  x=plot_res$mean_skew_obs,
  mean=plot_res$mean_bias,
  lower=plot_res$ci_low_bias,
  upper=plot_res$ci_high_bias
)


## plot

pC<-ggplot()+
  
  #population
  geom_line(data=pop_line,aes(x=x,y=mean),
            linewidth=0.7,color="#D5B60A")+
  geom_point(data=pop_line,aes(x=x,y=mean),
             size=2,color="#D5B60A")+
  
  #true cohort bias
  geom_line(data=true_line,aes(x=x,y=mean),
            linewidth=0.7,color="purple4")+
  geom_errorbar(data=true_line,
                aes(x=x,ymin=lower,ymax=upper),
                width=0.2,linewidth=0.5,
                alpha=0.7,color="purple4")+
  geom_point(data=true_line,aes(x=x,y=mean),
             size=2,color="purple4")+
  
  #estimated bias
  geom_line(data=est_line,aes(x=x,y=mean),
            linewidth=0.7,color="black")+
  geom_errorbar(data=est_line,
                aes(x=x,ymin=lower,ymax=upper),
                width=0.2,linewidth=0.5,
                alpha=0.7,color="black")+
  geom_point(data=est_line,aes(x=x,y=mean),
             size=2,color="black")+
  
  theme_classic(base_size=10)+
  theme(axis.title=element_blank(),
        axis.text=element_text(size=10),
        legend.position="none")+
  xlim(1,10.25)+
  ylim(-0.1,1)


####################################################################################
### final plot

top_fig3<-plot_grid(pA,pB,pC,nrow=1)

ggsave("top_fig3.svg",plot=top_fig3,
       width=921,height=227,units="px",
       device="svg")