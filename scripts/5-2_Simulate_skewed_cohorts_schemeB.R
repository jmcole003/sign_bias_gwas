## Simulation of skewed cohorts (scheme B)

# n_pop: population size
# n_pairs: number of matched pairs (increasing/decreasing)
# h2_target: target SNP h2
# n_samp: cohort size
# n_rep: number of reps per target skew
# block_size: pseudo ld block size
# maf_max: MAF threshold 
# maf_min_pop: minimum population MAF
# maf_max_pop: maximum population MAF
# tail_frac: fraction of the population in each tail mode
# maf_shift: allele frequency shift among modes
# center_lambda: sampling concentration ( center mode)
# right_lambda: sampling concentration (right mode)
# band_sd_center: width of band in center (SD units)
# band_sd_right: width of band in right mode
# right_frac_map: fraction sampled from the right mode


library(ggplot2)
library(ashr)
library(cowplot)

# params
n_pop<-2000000
n_pairs<-2000
h2_target<-0.6
n_rep<-20
n_samp<-10000
maf_min_pop<-0.0001
maf_max_pop<-0.1
tail_frac<-0.1
#tail_frac<-0.2
maf_shift<-0.5
#maf_shift<-0.1
#maf_shift<-0.3
#maf_shift<-0.8
pop_seed<-42
base_seed<-1
gwas_chunk<-200
block_size<-6
ci_level<-0.95
maf_max<-0.01
center_lambda<-200
right_lambda<-2
band_sd_center<-0.1
band_sd_right<-0.037

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
build_population <- function(n_pop,n_pairs,
                           h2=0.6,maf_min=0.0001, 
                           maf_max=0.1,
                           tail_frac, maf_shift,
                           seed=42){
  
  set.seed(seed)
  
  # number of snps and mafs
  n_snps <- n_pairs*2
  p_pair<-10^(runif(n_pairs, 
                    log10(maf_min), 
                    log10(maf_max)))
  p<-rep(p_pair, each=2) 

  bmag <- sqrt(h2/sum(2*p*(1-p)))
              
  # equal num of pos and neg              
  beta <-c(rep(bmag, n_pairs), rep(-bmag, n_pairs)) 
  
  #assign pop modes
  n_left<-round(n_pop*tail_frac)
  n_right<-round(n_pop*tail_frac)
  
  mode<-rep("center",n_pop)
  
  left<-sample.int(n_pop,n_left)
  mode[left]<-"left"
  
  remaining<-which(mode=="center")
  right<-sample(remaining,n_right)
  mode[right]<-"right"
  mode<-factor(mode,levels=c("left","center","right"))
  
  #genotypes
  G_raw<-matrix(as.raw(0),nrow=n_pop,ncol=n_snps)
  #dim(G_raw)
  
  gval<-rep(0,n_pop)
  
  for(i in 1:n_snps){
    
    p0<-p[i]
    
    #shift allele frequency according to mode and effect direction
    if(beta[i]>0){
      p_left<-p0*(1-maf_shift)
      p_right<-p0*(1+maf_shift)
    } else {
      p_left<-p0*(1+maf_shift)
      p_right<-p0*(1-maf_shift)
    }
    
    p_ind<-ifelse(mode=="left",p_left,
                  ifelse(mode=="right",p_right,p0))
    
    G<-rbinom(n_pop,2,p_ind)
    
    G_raw[,i]<-as.raw(G)
    gval<-gval+beta[i]*G
  }
  
  #add noise
  eps<-rnorm(n_pop,sd=1-h2)
  Y_pop<-gval+eps
  
  #out
  out<-list("N_pop"=n_pop,
            "M"=n_snps,
            "p"=p,
            "beta"=beta,
            "G_raw"=G_raw,
            "Y_pop"=Y_pop,
            "gval"=gval,
            "mode"=mode,
            "h2"=h2,
            "bmag"=bmag,
            "tail_frac"=tail_frac,
            "maf_shift"=maf_shift)
  
  return(out)
}


##cohort sampler from center and right modes
sample_cohort<-function(pop,n_samp,right_frac=0.05,
                        center_lambda=200,right_lambda=2,
                        band_sd_center=0.1,band_sd_right=0.03,
                        seed=1){
  
  set.seed(seed)
  
  #number sampled from each mode
  n_right<-round(n_samp*right_frac)
  n_center<-n_samp-n_right
  
  
  #center mode
  center_inds<-which(pop$mode=="center")
  center_y<-pop$Y_pop[center_inds]
  
  center_mean<-median(center_y)
  center_sd<-sd(center_y)
  
  #only sample near center of mode
  center_keep<-abs(center_y-center_mean)<=band_sd_center*center_sd
  center_inds<-center_inds[center_keep]
  center_y<-center_y[center_keep]
  
  #sampling weights
  center_weights<-exp(-center_lambda*
                        ((center_y-center_mean)/center_sd)^2)
  
  
  #right mode
  right_inds<-which(pop$mode=="right")
  right_y<-pop$Y_pop[right_inds]
  
  right_mean<-median(right_y)
  right_sd<-sd(right_y)
  
  #only sample near center of mode
  right_keep<-abs(right_y-right_mean)<=band_sd_right*right_sd
  right_inds<-right_inds[right_keep]
  right_y<-right_y[right_keep]
  
  #sampling weights
  right_weights<-exp(-right_lambda*
                       ((right_y-right_mean)/right_sd)^2)
  
  
  #sample individuals
  center_sample<-sample(center_inds,n_center,
                        replace=FALSE,prob=center_weights)
  
  right_sample<-sample(right_inds,n_right,
                       replace=FALSE,prob=right_weights)
  
  idx<-c(center_sample,right_sample)
  Y_samp<-pop$Y_pop[idx]
  
  
  #out
  out<-list("idx"=idx,
            "Y_samp"=Y_samp,
            "skew"=skew3(Y_samp),
            "right_frac"=right_frac)
  
  return(out)
}

# GWAS (linreg)
gwas_all_snps<-function(pop,idx,chunk=200){
  
  n<-length(idx)
  
  #phenotype
  y<-pop$Y_pop[idx]
  y<-y-mean(y)
  Syy<-sum(y^2)
  
  n_snps<-ncol(pop$G_raw)
  
  #save results
  betahat<-rep(NA,n_snps)
  se<-rep(NA,n_snps)
  af<-rep(NA,n_snps)
  maf<-rep(NA,n_snps)
  mac<-rep(NA,n_snps)
  
  #run SNPs in chunks
  for(start_snp in seq(1,n_snps,by=chunk)){
    
    end_snp<-min(start_snp+chunk-1,n_snps)
    snps<-start_snp:end_snp
    
    #genotypes
    G<-matrix(as.integer(pop$G_raw[idx,snps,drop=FALSE]),
              nrow=n,ncol=length(snps))
    
    #allele counts/frequencies
    allele_count<-colSums(G)
    
    af_snp<-allele_count/(2*n)
    maf_snp<-pmin(af_snp,1-af_snp)
    mac_snp<-pmin(allele_count,2*n-allele_count)
    
    af[snps]<-af_snp
    maf[snps]<-maf_snp
    mac[snps]<-mac_snp
    
    #linear regression
    gmean<-allele_count/n
    Sxx<-colSums(G^2)-n*(gmean^2)
    Sxy<-as.numeric(crossprod(G,y))
    
    b<-Sxy/Sxx
    
    #effect for minor allele
    flip<-af_snp>0.5
    b[flip]<--b[flip]
    
    #standard error
    SSE<-Syy-(Sxy^2)/Sxx
    sigma2<-SSE/(n-2)
    s<-sqrt(sigma2/Sxx)
    
    #save
    betahat[snps]<-b
    se[snps]<-s
  }
  
  #z score
  z<-betahat/se
  
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

#fit ASH
run_ash_sb<-function(betahat,se){
  
  fit<-ashr::ash(betahat=betahat,sebetahat=se,
                 mixcompdist="normal",method="fdr")
  
  sb<-fit$result$PositiveProb-fit$result$NegativeProb
  
  return(as.numeric(sb))
}

#cohort sign bias
cohort_sign_bias_from_sb<-function(sb){
  bias<-sum(sb)/sum(abs(sb))
  return(bias)
}

#true effect sign, oriented to minor allele
true_minor_sign<-function(af,beta){
  true_sign<-sign(beta)
  true_sign[af>0.5]<--true_sign[af>0.5]
  return(true_sign)
}

#lowest p-value SNP from each pseudo-block
min_p_per_block<-function(gwas,inds,block_size=6){
  
  z<-gwas$z[inds]
  p<-2*pnorm(-abs(z))
  
  block<-((gwas$snp[inds]-1)%/%block_size)+1
  blocks<-sort(unique(block))
  
  winners<-rep(NA,length(blocks))
  
  for(i in 1:length(blocks)){
    in_block<-which(block==blocks[i])
    min_p<-which.min(p[in_block])
    winners[i]<-inds[in_block[min_p]]
  }
  
  return(winners)
}

#fit ASH and calculate sign bias for one cohort
one_cohort_bias<-function(pop,samp,gwas_chunk=200,
                          block_size=6,maf_max=Inf){
  
  #GWAS
  gwas<-gwas_all_snps(pop,idx=samp$idx,chunk=gwas_chunk)
  
  #SNPs that can be used for ASH
  use_ash<-!is.na(gwas$betahat) & !is.na(gwas$se) &
    gwas$se>0 & gwas$mac>=1
  
  sb<-rep(NA,nrow(gwas))
  sb[use_ash]<-run_ash_sb(gwas$betahat[use_ash],gwas$se[use_ash])
  
  #rare SNPs used for sign bias
  use<-use_ash & gwas$maf>0 & gwas$maf<=maf_max &
    !is.na(sb) & !is.na(gwas$z)
  
  inds<-which(use)
  
  #true sign bias
  beta_true<-pop$beta[gwas$snp[inds]]
  true_sign<-true_minor_sign(gwas$af[inds],beta_true)
  true_bias<-mean(true_sign)
  
  #lowest p SNP from each pseudo-block
  winners<-min_p_per_block(gwas,inds,block_size=block_size)
  
  #estimated sign bias
  bias_sig<-cohort_sign_bias_from_sb(sb[winners])
  
  #out
  out<-data.frame(skew_obs=samp$skew,
                  bias_sig=bias_sig,
                  true_bias_sample=true_bias,
                  n_used=length(winners))
  
  return(out)
}

####################################################################################
##run simulations for each target skew

right_frac_map <-c(
  `0`= 0.00,
  `1`= 0.276,
  `2`=0.146,
  `3`= 0.0835,
  `4`= 0.0523,
  `5`= 0.0355,
  `6`= 0.0254,
  `7`= 0.0189,
  `8`= 0.0146,
  `9`= 0.01165,
  `10` = 0.0094
)

simulate_signbias_vs_skew<-function(pop,right_frac_map,n_rep=20,n_samp=10000,
                                    center_lambda=200,right_lambda=2,
                                    band_sd_center=0.1,band_sd_right=0.03,
                                    gwas_chunk=200,block_size=6,
                                    maf_max=Inf,ci_level=0.95,
                                    base_seed=1){
  
  targets<-as.numeric(names(right_frac_map))
  
  #### run replicates
  
  reps_list<-list()
  k<-1
  
  for(i in 1:length(targets)){
    
    target<-targets[i]
    right_frac<-right_frac_map[i]
    
    for(r in 1:n_rep){
      
      seed<-base_seed+100000*target+r
      
      samp<-sample_cohort(
        pop,
        n_samp=n_samp,
        right_frac=right_frac,
        center_lambda=center_lambda,
        right_lambda=right_lambda,
        band_sd_center=band_sd_center,
        band_sd_right=band_sd_right,
        seed=seed
      )
      
      res<-one_cohort_bias(
        pop,samp,
        gwas_chunk=gwas_chunk,
        block_size=block_size,
        maf_max=maf_max
      )
      
      reps_list[[k]]<-data.frame(
        target_sk=target,
        rep=r,
        right_frac=right_frac,
        skew_obs=res$skew_obs,
        bias_sig=res$bias_sig,
        true_bias_sample=res$true_bias_sample,
        n_used=res$n_used,
        seed=seed
      )
      
      k<-k+1
    }
  }
  
  reps<-do.call(rbind,reps_list)
  
  #### summarize
  
  alpha<-1-ci_level
  sum_list<-list()
  
  for(i in 1:length(targets)){
    
    target<-targets[i]
    x<-reps[reps$target_sk==target,]
    
    n_ok<-sum(!is.na(x$bias_sig))
    
    mean_bias<-mean(x$bias_sig,na.rm=TRUE)
    sd_bias<-sd(x$bias_sig,na.rm=TRUE)
    se_bias<-sd_bias/sqrt(n_ok)
    
    mean_true<-mean(x$true_bias_sample,na.rm=TRUE)
    sd_true<-sd(x$true_bias_sample,na.rm=TRUE)
    se_true<-sd_true/sqrt(n_ok)
    
    tcrit<-qt(1-alpha/2,df=n_ok-1)
    
    sum_list[[i]]<-data.frame(
      target_sk=target,
      right_frac=right_frac_map[i],
      n_rep=n_ok,
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
  
  summ<-do.call(rbind,sum_list)
  
  #out
  out<-list("reps"=reps,
            "summary"=summ)
  
  return(out)
}

####################################################################################
## RUN

#simulate population
pop<-build_population(
  n_pop=n_pop,
  n_pairs=n_pairs,
  h2=h2_target,
  maf_min=maf_min_pop,
  maf_max=maf_max_pop,
  tail_frac=tail_frac,
  maf_shift=maf_shift,
  seed=pop_seed
)

#run simulations
out<-simulate_signbias_vs_skew(
  pop,
  right_frac_map=right_frac_map,
  n_rep=n_rep,
  n_samp=n_samp,
  center_lambda=center_lambda,
  right_lambda=right_lambda,
  band_sd_center=band_sd_center,
  band_sd_right=band_sd_right,
  gwas_chunk=gwas_chunk,
  block_size=block_size,
  maf_max=maf_max,
  ci_level=ci_level,
  base_seed=base_seed
)

#example cohort for plotting
example_cohort<-sample_cohort(
  pop,
  n_samp=n_samp,
  right_frac=0.0094,
  center_lambda=center_lambda,
  right_lambda=right_lambda,
  band_sd_center=band_sd_center,
  band_sd_right=band_sd_right,
  seed=1
)

####################################################################################
##### Plotting

## Population distribution

dens_pop<-density(pop$Y_pop)
pop_plot<-data.frame(x=dens_pop$x,y=dens_pop$y)

mode_centers<-data.frame(
  mode=levels(pop$mode),
  x=as.numeric(tapply(pop$Y_pop,pop$mode,mean))
)

pD<-ggplot(pop_plot,aes(x,y))+
  geom_area(fill="#FDE992",color=NA)+
  geom_line(color="black",linewidth=0.6)+
  geom_vline(data=mode_centers,aes(xintercept=x),
             linewidth=0.8,linetype=2)+
  theme_classic()+
  theme(axis.title=element_blank(),
        axis.text=element_text(size=10))

## Example cohort distribution

y_example<-example_cohort$Y_samp
dens_cohort<-density(y_example,adjust=3)
cohort_plot<-data.frame(x=dens_cohort$x,y=dens_cohort$y)

pE<-ggplot(cohort_plot,aes(x,y))+
  geom_area(fill="#EFEFEF",color=NA)+
  geom_line(color="black",linewidth=0.6)+
  theme_classic()+
  theme(axis.title=element_blank(),
        axis.text=element_text(size=10))

####################################################################################
## Sign bias

res_sim<-out$summary

#data for each line
pop_line<-data.frame(
  x=res_sim$mean_skew_obs,
  mean=0
)

true_line<-data.frame(
  x=res_sim$mean_skew_obs,
  mean=res_sim$mean_true,
  lower=res_sim$ci_low_true,
  upper=res_sim$ci_high_true
)

est_line<-data.frame(
  x=res_sim$mean_skew_obs,
  mean=res_sim$mean_bias,
  lower=res_sim$ci_low_bias,
  upper=res_sim$ci_high_bias
)

pF<-ggplot()+
  
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
                width=0.15,linewidth=0.5,color="purple4")+
  geom_point(data=true_line,aes(x=x,y=mean),
             size=2,color="purple4")+
  
  #estimated sign bias
  geom_line(data=est_line,aes(x=x,y=mean),
            linewidth=0.7,color="black")+
  geom_errorbar(data=est_line,
                aes(x=x,ymin=lower,ymax=upper),
                width=0.15,linewidth=0.5,color="black")+
  geom_point(data=est_line,aes(x=x,y=mean),
             size=2,color="black")+
  
  theme_classic(base_size=10)+
  theme(axis.title=element_blank(),
        axis.text=element_text(size=10),
        legend.position="none")+
  xlim(1,10.1)+
  ylim(-0.05,0.4)

####################################################################################
### final plot

bottom_fig3<-plot_grid(pD,pE,pF,nrow=1)

ggsave("bottom_fig3.svg",plot=bottom_fig3,
       width=921,height=227,units="px",
       device="svg")