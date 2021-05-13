# This is the main function for our proposed method in Cox mediation analysis
#' High-dimensional mediation analysis for survival data
#' 
#' \code{survHIMA} is used to estimate and test high-dimensional mediation effects for survival data.
#' 
#' @param X a vector of exposure. 
#' @param Z a matrix of adjusting covariates. Rows represent samples, columns represent variables.
#' @param M a \code{data.frame} or \code{matrix} of high-dimensional mediators. Rows represent samples, columns 
#' represent variables.
#' @param y a vector of observed failure times.
#' @param status a vector of censoring indicator (\code{status = 1}: uncensored; \code{status = 0}: censored)
#' @param FDRcut FDR cutoff applied to define and select significant mediators. Default = \code{0.05}. 
#' 
#' @return A data.frame containing mediation testing results of selected mediators (FDR <\code{0.05}). 
#' \itemize{
#'     \item{ID: }{index of selected significant mediator.}
#'     \item{alpha: }{coefficient estimates of exposure (X) --> mediators (M).}
#'     \item{alpha_se: }{standard error for alpha.}
#'     \item{beta: }{coefficient estimates of mediators (M) --> outcome (Y) (adjusted for exposure).}
#'     \item{beta_se: }{standard error for beta}
#'     \item{p.joint: }{joint p-value of selected significant mediator.}
#' }
#'
#' @examples
#' ## Generate simulated survival data
#' n <- 300  # sample size
#' c_unif <- 150 # CR= 0.20
#' 
#' ##
#' p <- 100 # the dimension of mediators
#' q <- 2   # the dimension of covariates
#' ##
#' 
#' sigma_e <- matrix(0.25, p, p)
#' diag(sigma_e) <- 1
#' sigma_e[1, 3] <- 0.8
#' sigma_e[3, 1] <- sigma_e[1, 3]
#' sigma_e[2, 4] <- 0.8
#' sigma_e[4, 2] <- 0.8
#' 
#' ##
#' beta <- matrix(0, 1, p)
#' beta[1]  <- 0.6
#' beta[3]  <- 0.8
#' beta[5]  <- 0.55
#' 
#' ##
#' alpha <- matrix(0, 1, p)
#' alpha[1]  <- 0.6
#' alpha[3]  <- 0.8
#' alpha[5]  <- 0.8
#' 
#' ##
#' gamma <- matrix(0.5, 1, q)
#' eta <- matrix(0.3, p, q)
#' r <- matrix(0.5, 1, 1)
#' 
#' ##
#' X <- matrix(rnorm(n, mean = 0, sd = 2), n, 1) # expoure
#' Z <- matrix(rnorm(n * q, mean = 0, sd = 2), n, q) # covariates
#' mu <- matrix(0, p, 1)
#' e <- MASS::mvrnorm(n, mu, sigma_e)  # the error terms
#' 
#' M <- matrix(0, n, p)
#' M <- X%*%(alpha) + Z%*%t(eta) + e
#' MZ <- cbind(M, Z, X)
#' 
#' beta_gamma <- cbind(beta,gamma,r)
#' 
#' ## generate the filure time T 
#' u <- runif(n, 0, 1)
#' T <- matrix(0, n, 1) 
#' for (i in 1:n){
#'   T[i] <- -log(1 - u[i])*exp(-sum(beta_gamma*MZ[i,]))
#' }
#'
#' ## generate censoring time 0.45 censoring rate
#' C <- runif(n, min = 0, max = c_unif)  
#' 
#' status <- as.integer(T < C)
#' 
#' ## the observed failure time
#' y <- apply(cbind(C,T), 1, min) 
#' 
#' fit <- survHIMA(X, Z, M, y, status)
#' 
#' ## ID = 1, 3, 5, are selected (FDR <0.05)
#' fit
#' 
#' @export
survHIMA <- function(X, Z, M, y, status, FDRcut = 0.05){

  MZ <- cbind(M,Z,X)
  n <- length(X)
  p <- dim(M)[2]
  q <- dim(Z)[2]
  
  #########################################################################
  ################################ STEP 1 #################################
  #########################################################################
  message("Step 1: Mediators screening ...", "     (", Sys.time(), ")")
  
  d_0 <- 1*round(n/log(n))
  beta_SIS <- matrix(0,1,p)
  for (i in 1:p){
    ID_S <- c(i, (p+1):(p+q+1))
    MZ_SIS <- MZ[,ID_S]
    fit <- survival::coxph(survival::Surv(y, status) ~ MZ_SIS)
    beta_SIS[i] <- fit$coefficients[1]
  }

  alpha_SIS <- matrix(0,1,p)
  XZ <- cbind(X,Z)
  for (i in 1:p){
    fit_a  <- lsfit(XZ,M[,i],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    alpha_SIS[i] <- est_a
  }
  
  ab_SIS <- alpha_SIS*beta_SIS
  ID_SIS  <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[d_0])

  d <- length(ID_SIS)
  
  #########################################################################
  ################################ STEP 2 #################################
  #########################################################################
  message("Step 2: De-biased Lasso estimates ...", "     (", Sys.time(), ")")
  
  P_beta_SIS <- matrix(0,1,d)
  beta_DLASSO_SIS_est <- matrix(0,1,d)
  beta_DLASSO_SIS_SE <- matrix(0,1,d)
  MZ_SIS <- MZ[,c(ID_SIS, (p+1):(p+q+1))]
  MZ_SIS_1 <- t(t(MZ_SIS[,1]))
  
  for (i in 1:d){
    V <- MZ_SIS
    V[,1] <- V[,i]
    V[,i] <- MZ_SIS_1
    LDPE_res <- LDPE_func(ID = 1, X = V, time = y, status = status)
    beta_LDPE_est <- LDPE_res[1]
    beta_LDPE_SE  <- LDPE_res[2]
    V1_P <- abs(beta_LDPE_est)/beta_LDPE_SE
    P_beta_SIS[i] <- 2*(1-pnorm(V1_P,0,1))
    beta_DLASSO_SIS_est[i] <- beta_LDPE_est
    beta_DLASSO_SIS_SE[i] <- beta_LDPE_SE
  }
  
  ## the following is the code for the estimation of alpha
  alpha_SIS_est <- matrix(0,1,d)
  alpha_SIS_SE <- matrix(0,1,d)
  P_alpha_SIS <- matrix(0,1,d)
  XZ <- cbind(X,Z)
  
  for (i in 1:d){
    fit_a  <- lsfit(XZ,M[,ID_SIS[i]],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    se_a <- ls.diag(fit_a)$std.err[2]
    sd_1 <- abs(est_a)/se_a
    P_alpha_SIS[i] <- 2*(1-pnorm(sd_1,0,1))  ## the SIS for alpha
    alpha_SIS_est[i] <- est_a
    alpha_SIS_SE[i] <- se_a
  }
  
  #########################################################################
  ################################ STEP 3 #################################
  #########################################################################
  message("Step 3: Multiple-testing procedure ...", "     (", Sys.time(), ")")
  
  PA <- cbind(t(P_alpha_SIS), t(P_beta_SIS))
  P_value <- apply(PA,1,max)  # the joint p-values for SIS variable

  ## the multiple-testing  procedure
  N0 <- dim(PA)[1]*dim(PA)[2]
  
  input_pvalues <- PA + matrix(runif(N0,0,10^{-10}),dim(PA)[1],2)
  nullprop <- HDMT::null_estimation(input_pvalues,lambda=0.5)
  fdrcut  <- HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10, nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
  
  ID_fdr <- which(fdrcut <= FDRcut)

  if (length(ID_fdr) > 0){
    alpha_hat <- alpha_SIS_est[ID_fdr]
    alpha_est <- alpha_SIS_SE[ID_fdr]
    beta_hat <- beta_DLASSO_SIS_est[ID_fdr]
    beta_est <- beta_DLASSO_SIS_SE[ID_fdr]
    ID <- ID_SIS[ID_fdr]
    P_max <- P_value[ID_fdr]
  }
  
  out_result <- data.frame(ID = ID, 
                           alpha = alpha_hat, 
                           alpha_se = alpha_est, 
                           beta = beta_hat, 
                           beta_set = beta_est,
                           pvalue = P_max)
  
  return(out_result)
}
