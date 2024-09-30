#' @title Main function for heritability adjustment
#' @description heritability adjustment
#' @param path working path
#' @param mean_shift mean shift of the phenotype of interest, between the sample of participant (UKBB) and the population, standardized in the sample of participants ((mean_participants-mean_population)/SE_in_participants)
#' @param trait_name Name of the phenotype
#' @param trait_binary T/F, T for binary traits, F for the other trait types
#' @param K disease prevalence. Need to be specified for 'binary' phenotypes to derive observed-scale heritability
#' @param h2x heritability of participation liability score. Default=0.125 in UKBB
#' @param alpha participation rate. Default=0.055 in UKBB
#' @param lambda recurrence ratio. Default=2 in UKBB
#' @param shift_sd standard deviations of the mean shift, default=1/sqrt(HSE_sample_size)
#' @import data.table stats utils mvtnorm readr
#' @export


h2_PB_adjust <- function(path, mean_shift, trait_name='trait1',
                         trait_binary=F, K=1, h2x=0.125, alpha=0.055,
                         lambda=2, shift_sd=1/sqrt(20208)){
  set.seed(666)
  #path <- result_path
  system(paste0("grep ' ", path,
  "sumstats/PB.sumstats.gz  ", path, "sumstats/", trait_name, ".sumstats.gz' ",path,'results_',trait_name,"/res_rg.log  > ",path,'/results_',trait_name,"/res_rg_tab.log"))
  a=read_log(paste0(path,'results_',trait_name,'/res_rg_tab.log'), progress=F)
  if(a$X7=='NA'){
    return(NA)
  }
  h22_est <-  a$X7
  if(h22_est<0){
    return(NA)
  }
  h22_est.se <- a$X8
  if(trait_binary){
    P <- K
    h22_est <- cvt.h2(K,P,h22_est)
    h22_est.se <- cvt.h2(K,P,h22_est.se)
  }

  h22_est_pv <- pnorm(-abs(a$X7/a$X8))
  #h2.se <- c(h2.se, a$X8)

  gcor_est <- a$X3
  #h21_est <- 0.125

  correct0 <- get_correct_diff_h2_jack_noldsc(h21=h2x,h22_est=h22_est,rg_est=gcor_est, shift=mean_shift)
  nom <- data.frame(fread(paste0(path,'/results_',trait_name,'/nomvalues.txt'))) ### gcov
  #denom <- fread('denomvalues.txt') ### sqrt(h21*h22)
  hsq1 <- data.frame(fread(paste0(path,'/results_',trait_name,'/hsq1.txt')))
  hsq2 <- data.frame(fread(paste0(path,'/results_',trait_name,'/hsq2.txt')))
  hsq1.obs <- hsq1/ get_delta_prime_obs(alpha=alpha,lam=lambda,alpha_3=1/3)^2
  rg_est <- nom/sqrt(hsq1*hsq2)
  shift.sd <- rnorm(length(hsq1),mean=mean_shift,sd=shift_sd)
  correct.jack <- get_correct_diff_h2_jack_noldsc(hsq1.obs,hsq2,rg_est,shift=shift.sd)
  nblock <- nrow(nom)
  temp.h2 <- rep(nblock*unlist(correct0$h22_correct),nblock)-(nblock-1)*unlist(correct.jack$h22_correct)
  h22.se <- sqrt(var(temp.h2)/nblock)
  temp.rg <- rep(nblock*unlist(correct0$rg_correct),nblock)-(nblock-1)*unlist(correct.jack$rg_correct)
  rg.se <-  sqrt(var(temp.rg)/nblock)
  temp.re <- rep(nblock*unlist(correct0$re_correct),nblock)-(nblock-1)*unlist(correct.jack$re_correct)
  re.se <-  sqrt(var(temp.re)/nblock)

  return(list(h2y_adj=correct0$h22_correct,h2y_se_adj=h22.se,
              rhog_adj=correct0$rg_correct,
              rhog_se_adj=rg.se,
              rhoe_adj=correct0$re_correct,
              rhoe_se_adj=re.se,
              h2y_est=h22_est,
              h2y_se_est=h22_est.se,
              rhog_est=a$X3,
              rhog_se_est=a$X4))

}
#h2_PB_adjust(path,mean_shift,trait_name)


