#' @title Main function for the adjustment of genetic correlation between 2 phenotypes
#' @description heritability adjustment
#' @param path working path
#' @param mean_shift1 mean shift of the phenotype 1 of interest, between the sample of participant (UKBB) and the population, standardized in the sample of participants ((mean_participants-mean_population)/SE_in_participants)
#' @param mean_shift2 mean shift of the phenotype 2 of interest, between the sample of participant (UKBB) and the population, standardized in the sample of participants ((mean_participants-mean_population)/SE_in_participants)
#' @param trait_name1 Name of the phenotype 1
#' @param trait_name2 Name of the phenotype 2
#' @param trait_binary1 trait type of phenotype 1, T/F, T for binary traits, F for the other trait types
#' @param trait_binary2 trait type of phenotype 2, T/F, T for binary traits, F for the other trait types
#' @param K1 disease prevalence for phenotype 1. Need to be specified for 'binary' phenotypes to derive observed-scale heritability
#' @param K2 disease prevalence for phenotype 2. Need to be specified for 'binary' phenotypes to derive observed-scale heritability
#' @param h2x heritability of participation liability score. Default=0.125 in UKBB
#' @param alpha participation rate. Default=0.055 in UKBB
#' @param lambda recurrence ratio. Default=2 in UKBB
#' @param shift_sd1 standard deviations of the mean shift of phenotype 1, default=1/sqrt(HSE_sample_size)
#' @param shift_sd2 standard deviations of the mean shift of phenotype 2, default=1/sqrt(HSE_sample_size)
#' @import data.table stats utils mvtnorm readr
#' @export

gcor_PB_adjust <- function(path, mean_shift1, mean_shift2,
                           trait_name1='trait1', trait_name2='trait1',
                           trait_binary1=F, trait_binary2=F,
                           K1=1, K2=1, h2x=0.125, alpha=0.055, lambda=2,
                           shift_sd1=1/sqrt(20208), shift_sd2=1/sqrt(20208)){


  set.seed(666)
  # trait <- trait2
  # shift <- shift2
  if(T){
    # system(paste0("grep ' ", path,
    #               "sumstats/PB.sumstats.gz  ", path, "sumstats/", trait_name1, ".sumstats.gz' ",path,'results_',trait_name1,"/res_rg.log  > ",path,'/results_',trait_name1,"/res_rg_tab.log"))
    system(paste0("grep -A 2 'Summary of Genetic Correlation Results' ",path,'results_',trait_name,"/res_rg.log  | tail -n 1 > ",path,'/results_',trait_name,"/res_rg_tab.log"))


    a=read_log(paste0(path,'results_',trait_name1,'/res_rg_tab.log'), progress=F)
    if(a$X7=='NA'){
      return(NA)
    }
    h22_est <-  a$X7
    if(h22_est<0){
      return(NA)
    }
    h22_est.se <- a$X8
    if(trait_binary1){
      P <- K <- K1
      h22_est <- cvt.h2(K,P,h22_est)
      h22_est.se <- cvt.h2(K,P,h22_est.se)
    }
    h22_est_pv <- pnorm(-abs(a$X7/a$X8))
    #h2.se <- c(h2.se, a$X8)
    gcor_est <- a$X3
    #h21_est <- 0.125
    correct02 <- get_correct_diff_h2_jack_noldsc(h21=0.125,h22_est=h22_est,rg_est=gcor_est, shift=mean_shift1)
  }
  h22_est00 <- h22_est
  # trait <- trait3
  # shift <- shift3
  if(T){
    system(paste0("grep ' ", path,
                  "sumstats/PB.sumstats.gz  ", path, "sumstats/", trait_name2, ".sumstats.gz' ",path,'results_',trait_name2,"/res_rg.log  > ",path,'/results_',trait_name2,"/res_rg_tab.log"))
    a=read_log(paste0(path,'results_',trait_name2,'/res_rg_tab.log'), progress=F)
    if(a$X7=='NA'){
      return(NA)
    }
    h22_est <-  a$X7
    if(h22_est<0){
      return(NA)
    }
    h22_est.se <- a$X8
    if(trait_binary2){
      P <- K <-  K2
      h22_est <- cvt.h2(K,P,h22_est)
      h22_est.se <- cvt.h2(K,P,h22_est.se)
    }

    h22_est_pv <- pnorm(-abs(a$X7/a$X8))
    #h2.se <- c(h2.se, a$X8)
    gcor_est <- a$X3
    #h21_est <- 0.125
    correct03 <- get_correct_diff_h2_jack_noldsc(h21=0.125,h22_est=h22_est,rg_est=gcor_est, shift=mean_shift2)
  }
  h23_est00 <- h22_est
  # pheno1 <- trait2
  # pheno2 <- trait3
  if(T){
    system(paste0("grep ' ", path,
                  "sumstats/", trait_name1,".sumstats.gz  ", path, "sumstats/", trait_name2, ".sumstats.gz' ",path,'results_',trait_name1, '_', trait_name2,"/res_rg.log  > ",path,'/results_',trait_name1,'_', trait_name2,"/res_rg_tab.log"))
    a=read_log(paste0(path,'results_',trait_name1, '_', trait_name2,'/res_rg_tab.log'), progress=F)
    phig_est <- a$X3
    phig_est_se <- a$X4
  }
  phiG_est00 <- phig_est*sqrt(h22_est00*h23_est00)
  phiG_correct00 <- get_correct_diff_gcor_jack_noldsc(0.125,
                                                      rho_T_correct2=correct02$rho_T_correct,
                                                      rho_T_correct3=correct03$rho_T_correct,
                                                      rho_G_correct2=correct02$rho_G_correct,
                                                      rho_G_correct3=correct03$rho_G_correct,
                                                      phiG_est=phiG_est00)
  phig_correct00 <- phiG_correct00/sqrt(correct02$h22_correct*correct03$h22_correct)

  #### jackknife
  # trait <- trait2
  # shift <- shift2
  if(T){
    nom2 <- data.frame(fread(paste0(path,'/results_',trait_name1,'/nomvalues.txt')) )### gcov
    #denom <- fread('denomvalues.txt') ### sqrt(h21*h22)
    hsq1 <- data.frame(fread(paste0(path,'/results_',trait_name1,'/hsq1.txt')))
    hsq2 <- data.frame(fread(paste0(path,'/results_',trait_name1,'/hsq2.txt')))
    hsq1.obs <- hsq1/ get_delta_prime_obs(alpha=alpha,lam=lambda,alpha_3=1/3)^2
  }
  # trait <- trait3
  # shift <- shift3
  if(T){
    #setwd(paste0('/home/songs/UKB_summary/20210211/clean.ldsc/summs/',trait))
    nom3 <- data.frame(fread(paste0(path,'/results_',trait_name2,'/nomvalues.txt'))) ### gcov
    #denom <- fread('denomvalues.txt') ### sqrt(h21*h22)
    hsq3 <- data.frame(fread(paste0(path,'/results_',trait_name2,'/hsq2.txt')))
  }
  # pheno1 <- trait2
  # pheno2 <- trait3
  if(T){
    #setwd(paste0('/home/songs/datasets/HSE/gcor/',pheno1,'_',pheno2))
    #nom23 <- fread('nomvalues.txt') ### gcov
    nom23 <- data.frame(fread(paste0(path,'/results_',trait_name1,'_',trait_name2, '/nomvalues.txt'))) ### gcov

  }
  rg_est2 <- nom2/sqrt(hsq1*hsq2)
  correct.jack2 <- get_correct_diff_h2_jack_noldsc(hsq1.obs,hsq2,rg_est2,shift=rnorm(length(hsq1),mean=mean_shift1,sd=shift_sd1))
  rg_est3 <- nom3/sqrt(hsq1*hsq3)
  correct.jack3 <- get_correct_diff_h2_jack_noldsc(hsq1.obs,hsq3,rg_est3,shift=rnorm(length(hsq1),mean=mean_shift2,sd=shift_sd2))

  phiG_correct.jack <- get_correct_diff_gcor_jack_noldsc(hsq1.obs,
                                                         rho_T_correct2=correct.jack2$rho_T_correct,
                                                         rho_T_correct3=correct.jack3$rho_T_correct,
                                                         rho_G_correct2=correct.jack2$rho_G_correct,
                                                         rho_G_correct3=correct.jack3$rho_G_correct,
                                                         phiG_est=nom23)

  phig_correct.jack <- phiG_correct.jack/sqrt(correct.jack2$h22_correct*correct.jack3$h22_correct)

  nblock <- nrow(nom23)
  temp.phiG <- rep(nblock*unlist(phiG_correct00),nblock)-
    (nblock-1)*unlist(phiG_correct.jack)
  phiG.se <- sqrt(var(temp.phiG)/nblock)
  temp.phig <- rep(nblock*unlist(phig_correct00),nblock)-
    (nblock-1)*unlist(phig_correct.jack)
  phig.se <- sqrt(var(temp.phig)/nblock)

  return(list(phig_adj=phig_correct00,
              phig_se_ajd=phig.se,
              phig_est=phig_est,
              phig_se_est=phig_est_se))


}





# gcor_PB_adjust(path='/n/home04/ssong0110/testRpackage/',
#                mean_shift1=0.438, mean_shift2=-0.138,
#                trait_name1='EA', trait_name2='BMI',
#                trait_binary1=F, trait_binary2=F,
#                K1=1, K2=1, h2x=0.125, alpha=0.055, lambda=2,
#                shift_sd1=1/sqrt(20208), shift_sd2=1/sqrt(20208))




