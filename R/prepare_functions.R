#' @import  stats utils mvtnorm

cvt.h2 <- function(K,P,h2){
  # K: prevalence
  # P: proportion of cases
  zv <- dnorm(qnorm(K))
  h2_liab <- h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
  return(h2_liab)
}

get_rho <- function(lam, alpha){
  rho <-  uniroot(opt, interval=c(-1,1), lam=lam, alpha=alpha)$root
  return(rho)
}
opt <- function(lam, alpha, rho){
  cut <- -qnorm(alpha,sd=1)
  ## lam * alpha^2 - P(sib1=1, sib2=1)
  ## find zero point
  lower <- rep(cut, 2)
  # upper <- rep(3, 5)
  corr <- diag(2)
  corr[lower.tri(corr)] <- rho
  corr[upper.tri(corr)] <- rho
  prob <- pmvnorm(lower=lower,
                  corr=corr)
  val <- lam*alpha^2-prob
  return(val)
}
get_lam1 <- function(cc,rho){
  upper <- c(-cc,-cc)
  # upper <- rep(3, 5)
  corr <- diag(2)
  #rho <- 0.5
  corr[lower.tri(corr)] <- rho
  corr[upper.tri(corr)] <- rho
  prob2 <- pmvnorm(upper=upper,
                   corr=corr)[1]
  ####
  lam1 <- dnorm(-cc)*pnorm(-cc*sqrt((1-rho)/(1+rho)))/prob2
  return(lam1)
}


get_delta_prime_obs <- function(alpha, lam, alpha_3=1/2){
  ## transformation:
  ## beta'=delta'*beta  (beta' is those derived with IBD)
  #delta <- approx_ratio(alpha,lam)
  delta <- get_ratio_shrink(alpha,lam) # f2-f0/fsam-f
  #tau_3 <- -qnorm(alpha_3)
  cc <- -qnorm(alpha)
  temp <- dnorm(-cc)/pnorm(-cc)
  temp_3 <-  sqrt(alpha_3*(1-alpha_3))
  delta_prime <- delta*temp*temp_3
  return(delta_prime)
}
get_ratio_shrink <- function(alpha, lam){
  cc <- -qnorm(alpha)
  rho <- get_rho(lam,alpha)
  #var_g <- seq(1e-4,0.1,by=0.001)
  # rhoU <- get_rhoU_new(lam,alpha,rho_g = 0.5, var_g=var_g)
  # rhoU
  return((get_lam1(cc,rho=rho))/
           (dnorm(-cc)/pnorm(-cc)))
}


get_correct_diff_gcor_jack_noldsc <- function(h21,rho_T_correct2,rho_T_correct3,
                                              rho_G_correct2,rho_G_correct3,phiG_est,
                                              alpha=0.055){
  if(T){
    cc <- -qnorm(alpha) ## cutoff
    ratio <- dnorm(cc)/alpha
    ksai <- ratio^2-cc*ratio
  }
  # h22_correct[k] <- h22_est[k]*(1-rho_T_true^2*ksai)+
  #   2*rho_T_true*rho_G_true*ksai-rho_T_true^2*h21*delta
  phiG_correct <- phiG_est*sqrt((1-rho_T_correct2^2*ksai)*(1-rho_T_correct3^2*ksai))-
    ksai^2*(rho_T_correct2*rho_T_correct3*h21)+
    ksai*(rho_T_correct3*rho_G_correct2+rho_T_correct2*rho_G_correct3)

  #phig_correct <- phiG_correct/sqrt(h22_correct*h23_correct)
  return(phiG_correct)
}




get_correct_diff_h2_jack_noldsc <- function(h21,h22_est,rg_est, alpha=0.055,shift){
  shift2 <- shift
  rg_est2 <- rg_est
  if(T){
    cc <- -qnorm(alpha) ## cutoff
    ratio <- dnorm(cc)/alpha
    ksai <- ratio^2-cc*ratio

    #rho2 <- list()
    rho_T_correct <- sign(shift2)/sqrt(ksai+dnorm(cc)^2/alpha^2/shift2^2)
    #rho_G_correct <- sqrt(1-rho_T_correct^2*ksai)*rg_est2*sqrt(h21*h22_est)  +ksai*rho_T_correct*h21
    rho_G_correct <- sqrt(1-rho_T_correct^2*ksai)*rg_est2*sqrt(h21*h22_est)*sqrt(1-ksai*h21)  +ksai*rho_T_correct*h21
    rho_E_correct <- rho_T_correct-rho_G_correct
  }
  # h22_correct[k] <- h22_est[k]*(1-rho_T_true^2*ksai)+
  #   2*rho_T_true*rho_G_true*ksai-rho_T_true^2*h21*delta
  h22_correct <- h22_est*(1-rho_T_correct^2*ksai)+
    2*rho_T_correct*rho_G_correct*ksai-rho_T_correct^2*h21*ksai^2-
    ksai/(1-ksai*h21)*(rho_G_correct-ksai*rho_T_correct*h21)^2
  rg_correct <- rho_G_correct/sqrt(h21* h22_correct)
  re_correct <- rho_E_correct/sqrt((1-h21)* (1-h22_correct))
  # rg_correct_new3[k] <- rho3$rho_G/sqrt(h21_est[k]* h23_correct_new[k])
  # phiG_correct <- phiG_est[k]*sqrt((1-rho2$rho_T^2*ksai)*(1-rho3$rho_T^2*ksai))-
  #   ksai^2*(rho2$rho_T*rho3$rho_T*h21_est[k])+
  #   ksai*(rho3$rho_T*rho2$rho_G+rho2$rho_T*rho3$rho_G)
  # phig_correct[k] <- phiG_correct/sqrt(h22_correct_new[k]*h23_correct_new[k])
  return(list(h22_correct=h22_correct,rg_correct=rg_correct,re_correct=re_correct,rho_T_correct=rho_T_correct,
              rho_G_correct=rho_G_correct))
}


