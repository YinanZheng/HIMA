# This is the main function for our proposed method for high-dimensional Cox mediation analysis
#' High-dimensional mediation analysis for survival data
#' 
#' \code{survHIMA} is used to estimate and test high-dimensional mediation effects for survival data.
#' 
#' @param X a vector of exposure. 
#' @param M a \code{data.frame} or \code{matrix} of high-dimensional mediators. Rows represent samples, columns 
#' represent mediator variables.
#' @param COV a matrix of adjusting covariates. Rows represent samples, columns represent variables. Can be \code{NULL}.
#' @param OT a vector of observed failure times.
#' @param status a vector of censoring indicator (\code{status = 1}: uncensored; \code{status = 0}: censored)
#' @param FDRcut HDMT pointwise FDR cutoff applied to select significant mediators. Default = \code{0.05}. 
#' @param scale logical. Should the function scale the data? Default = \code{TRUE}.
#' @param verbose logical. Should the function be verbose? Default = \code{FALSE}.
#' 
#' @return A data.frame containing mediation testing results of significant mediators (FDR <\code{FDRcut}). 
#' \itemize{
#'     \item{Index: }{mediation name of selected significant mediator.}
#'     \item{alpha_hat: }{coefficient estimates of exposure (X) --> mediators (M).}
#'     \item{alpha_se: }{standard error for alpha.}
#'     \item{beta_hat: }{coefficient estimates of mediators (M) --> outcome (Y) (adjusted for exposure).}
#'     \item{beta_se: }{standard error for beta.}
#'     \item{IDE: }{mediation (indirect) effect, i.e., alpha*beta.}
#'     \item{rimp: }{relative importance of the mediator.}
#'     \item{pmax: }{joint raw p-value of selected significant mediator (based on HDMT pointwise FDR method).}
#' }
#' 
#' @references Zhang H, Zheng Y, Hou L, Zheng C, Liu L. Mediation Analysis for Survival Data with High-Dimensional Mediators. 
#' Bioinformatics. 2021. DOI: 10.1093/bioinformatics/btab564. PMID: 34343267; PMCID: PMC8570823
#' 
#' @examples
#' \dontrun{
#' # Note: In the following example, M1, M2, and M3 are true mediators.
#' data(himaDat)
#' 
#' head(himaDat$Example3$PhenoData)
#' 
#' survHIMA.fit <- survHIMA(X = himaDat$Example3$PhenoData$Treatment,
#'                 M = himaDat$Example3$Mediator, 
#'                 COV = himaDat$Example3$PhenoData[, c("Sex", "Age")], 
#'                 OT = himaDat$Example3$PhenoData$Time, 
#'                 status = himaDat$Example3$PhenoData$Status, 
#'                 FDRcut = 0.05,
#'                 scale = FALSE, 
#'                 verbose = TRUE)
#' survHIMA.fit
#' }
#' 
#' @export
survHIMA <- function(X, M, COV = NULL, OT, status, 
                     FDRcut = 0.05, 
                     scale = TRUE, 
                     verbose = FALSE)
{
  
  X <- matrix(X, ncol = 1)
  M <- as.matrix(M)
  
  M_ID_name <- colnames(M)
  if(is.null(M_ID_name)) M_ID_name <- seq_len(ncol(M))
  
  n <- nrow(M)
  p <- nrow(M)
  
  if(is.null(COV))
  {q <- 0; MZ <- cbind(M,X)}
  else
  {COV <- as.matrix(COV); q <- dim(COV)[2]; MZ <- cbind(M,COV,X)}
  
  if(scale) 
  {
    MZ <- scale(MZ)
    if(verbose) message("Data scaling is completed.")
  }
  
  #########################################################################
  ################################ STEP 1 #################################
  #########################################################################
  message("Step 1: Sure Independent Screening ...", "     (", format(Sys.time(), "%X"), ")")
  
  d_0 <- 1*round(n/log(n))
  beta_SIS <- matrix(0,1,p)
  for (i in 1:p){
    ID_S <- c(i, (p+1):(p+q+1))
    MZ_SIS <- MZ[,ID_S]
    fit <- survival::coxph(survival::Surv(OT, status) ~ MZ_SIS)
    beta_SIS[i] <- fit$coefficients[1]
  }
  
  alpha_SIS <- matrix(0,1,p)
  XZ <- cbind(X,COV)
  for (i in 1:p){
    fit_a  <- lsfit(XZ,M[,i],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    alpha_SIS[i] <- est_a
  }
  
  ab_SIS <- alpha_SIS*beta_SIS
  ID_SIS  <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[min(p, d_0)])
  
  d <- length(ID_SIS)
  
  if(verbose) message("        Top ", d, " mediators are selected: ", paste0(M_ID_name[ID_SIS], collapse = ", "))
  
  #########################################################################
  ################################ STEP 2 #################################
  #########################################################################
  message("Step 2: De-biased Lasso estimates ...", "     (", format(Sys.time(), "%X"), ")")
  
  if(verbose)
  {
    if(is.null(COV)) 
    {message("        No covariate was adjusted.")} 
    else
    {message("        Adjusting for covariate(s): ", paste0(colnames(COV), collapse = ", "))}
  }
  
  ## estimation of beta
  P_beta_SIS <- matrix(0,1,d)
  beta_DLASSO_SIS_est <- matrix(0,1,d)
  beta_DLASSO_SIS_SE <- matrix(0,1,d)
  MZ_SIS <- MZ[,c(ID_SIS, (p+1):(p+q+1))]
  MZ_SIS_1 <- t(t(MZ_SIS[,1]))
  
  for (i in 1:d){
    V <- MZ_SIS
    V[,1] <- V[,i]
    V[,i] <- MZ_SIS_1
    LDPE_res <- LDPE_func(ID = 1, X = V, OT = OT, status = status)
    beta_LDPE_est <- LDPE_res[1]
    beta_LDPE_SE  <- LDPE_res[2]
    V1_P <- abs(beta_LDPE_est)/beta_LDPE_SE
    P_beta_SIS[i] <- 2*(1-pnorm(V1_P,0,1))
    beta_DLASSO_SIS_est[i] <- beta_LDPE_est
    beta_DLASSO_SIS_SE[i] <- beta_LDPE_SE
  }
  
  ## estimation of alpha
  alpha_SIS_est <- matrix(0,1,d)
  alpha_SIS_SE <- matrix(0,1,d)
  P_alpha_SIS <- matrix(0,1,d)
  XZ <- cbind(X,COV)
  
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
  message("Step 3: Multiple-testing procedure ...", "     (", format(Sys.time(), "%X"), ")")
  
  PA <- cbind(t(P_alpha_SIS), t(P_beta_SIS))
  P_value <- apply(PA, 1, max)  # the joint p-values for SIS variable
  
  ## the multiple-testing  procedure
  N0 <- dim(PA)[1]*dim(PA)[2]
  
  input_pvalues <- PA + matrix(runif(N0,0,10^{-10}),dim(PA)[1],2)
  nullprop <- null_estimation(input_pvalues,lambda=0.5)
  fdrcut  <- HDMT::fdr_est(nullprop$alpha00,
                           nullprop$alpha01,
                           nullprop$alpha10, 
                           nullprop$alpha1,
                           nullprop$alpha2,
                           input_pvalues,
                           exact=0)
  
  ID_fdr <- which(fdrcut <= FDRcut)
  
  IDE <- alpha_SIS_est[ID_fdr] * beta_DLASSO_SIS_est[ID_fdr]
  
  if (length(ID_fdr) > 0){
    out_result <- data.frame(Index = M_ID_name[ID_fdr], 
                             alpha_hat = alpha_SIS_est[ID_fdr], 
                             alpha_se = alpha_SIS_SE[ID_fdr], 
                             beta_hat = beta_DLASSO_SIS_est[ID_fdr], 
                             beta_se = beta_DLASSO_SIS_SE[ID_fdr],
                             IDE = IDE, 
                             rimp = abs(IDE)/sum(abs(IDE)) * 100, 
                             pmax = P_value[ID_fdr])
    if(verbose) message(paste0("        ", length(ID_fdr), " significant mediator(s) identified."))
  } else {
    if(verbose) message("        No significant mediator identified.")
    out_result = NULL
  }
  
  message("Done!", "     (", format(Sys.time(), "%X"), ")")
  
  return(out_result)
}
