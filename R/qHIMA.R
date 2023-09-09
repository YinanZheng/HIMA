# This is the main function for our proposed method for high-dimensional quantile mediation analysis
#' High-dimensional quantile mediation analysis
#' 
#' \code{qHIMA} is used to estimate and test high-dimensional quantile mediation effects.
#' 
#' @param X a vector of exposure. 
#' @param M a \code{data.frame} or \code{matrix} of high-dimensional mediators. Rows represent samples, columns 
#' represent mediator variables.
#' @param Y a vector of continuous outcome. Do not use data.frame or matrix.
#' @param Z a matrix of adjusting covariates. Rows represent samples, columns represent variables. Can be \code{NULL}.
#' @param penalty the penalty to be applied to the model (a parameter passed to function \code{conquer.cv.reg} in package \code{\link{conquer}}. 
#' Either 'MCP' (the default), 'SCAD', or 'lasso'.
#' @param topN an integer specifying the number of top markers from sure independent screening. 
#' Default = \code{NULL}. If \code{NULL}, \code{topN} will be \code{2*ceiling(n/log(n))}, where \code{n} is the sample size.
#' If the sample size is greater than topN (pre-specified or calculated), all mediators will be included in the test (i.e. low-dimensional scenario).
#' @param tau quantile level of outcome. Default = 0.5. A vector of tau is accepted.
#' @param scale logical. Should the function scale the data? Default = \code{TRUE}.
#' @param Bonfcut Bonferroni-corrected p value cutoff applied to define and select significant mediators. Default = \code{0.05}. 
#' @param verbose logical. Should the function be verbose? Default = \code{FALSE}.
#' @param ... other arguments.
#' 
#' @return A data.frame containing mediation testing results of selected mediators (Bonferroni-adjusted p value <\code{Bonfcut}). 
#' \itemize{
#'     \item{ID: }{index of selected significant mediator.}
#'     \item{alpha: }{coefficient estimates of exposure (X) --> mediators (M).}
#'     \item{alpha_se: }{standard error for alpha.}
#'     \item{beta: }{coefficient estimates of mediators (M) --> outcome (Y) (adjusted for exposure).}
#'     \item{beta_se: }{standard error for beta}
#'     \item{Bonferroni.p: }{statistical significance of the mediator (Bonferroni-corrected p value).}
#' }
#' 
#' @references Zhang H, Hong X, Zheng Y, Hou L, Zheng C, Wang X, Liu L. High-Dimensional Quantile Mediation Analysis with Application to a Birth 
#' Cohort Study of Mother–Newborn Pairs. Bioinformatics. 2023. (In press)
#' 
#' @examples
#' \dontrun{
#' data(Example5)
#' head(Example5$PhenoData)
#' 
#' qHIMA.fit <- qHIMA(X = Example5$PhenoData$Treatment,
#'                 M = Example5$Mediator, 
#'                 Y = Example5$PhenoData$Outcome, 
#'                 Z = Example5$PhenoData[, c("Sex", "Age")], 
#'                 Bonfcut = 0.05,
#'                 tau = c(0.3, 0.5, 0.7),
#'                 scale = FALSE, 
#'                 verbose = TRUE)
#' qHIMA.fit
#' }
#' 
#' @export
qHIMA <- function(X, M, Y, Z, penalty = c('MCP', "SCAD", "lasso"), 
                  topN = NULL, tau = 0.5, scale = TRUE, Bonfcut = 0.05, verbose = FALSE, ...){
  
  penalty <- match.arg(penalty)
  
  n <- nrow(M)
  p <- ncol(M)
  
  if(scale)
  {
    X <- scale(X)
    M <- scale(M)
    Z <- scale(Z)
  } else {
    X <- as.matrix(X)
    M <- as.matrix(M)
    Z <- as.matrix(Z)
  }
  
  if(is.null(topN)) {
    d <- ceiling(2 * n/log(n)) 
  } else {
    d <- topN  # the number of top mediators that associated with exposure (X)
  }
  
  d <- min(p, d) # if d > p select all mediators
  
  M_ID_name <- colnames(M)
  if(is.null(M_ID_name)) M_ID_name <- seq_len(p)
  
  #------------- Step 1: Mediator screening ---------------------------
  message("Step 1: Sure Independent Screening ...", "     (", format(Sys.time(), "%X"), ")")
  
  alpha_est <- matrix(0,1,p) # the OLS estimator of alpha
  alpha_SE  <- matrix(0,1,p) # the SE of alpha-OLS
  beta_SIS_est <- matrix(0,1,p) # the screening based estimator of beta
  beta_SIS_SE <- matrix(0,1,p) # the SE of beta_SIS_est
  
  XZ <- cbind(X,Z)
  for (k in 1:p){
    MXZ_k <- cbind(cbind(M[,k], X), Z)
    fit_M <- lsfit(XZ,M[,k],intercept = TRUE) # screening in the path x-M
    alpha_est[k] <- matrix(coef(fit_M))[2]
    alpha_SE[k] <- ls.diag(fit_M)$std.err[2]
  }
  
  T_sobel <- alpha_est
  ID_SIS <- which(-abs(T_sobel) <= sort(-abs(T_sobel))[d]) # the index set in Step 1 after the screening

  alpha_est_hat <- alpha_est
  
  if(verbose) message("    Top ", length(ID_SIS), " mediators are selected: ", paste0(M_ID_name[ID_SIS], collapse = ", "))
  
  
  #----------- Step 2: Penalized estimate in Quantile Regression model
  message("Step 2: Penalized estimate (", penalty, ") ...", "     (", format(Sys.time(), "%X"), ")")
  
  conf.names <- colnames(Z)
  if(verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
  
  MXZ_ID_SIS <- cbind(M[,ID_SIS],XZ)
  
  out_result <- NULL
  
  for (tau_temp in tau)
  {
    message("    Running penalized quantile regression with tau = ", tau_temp, " ...", "     (", format(Sys.time(), "%X"), ")")
    
    fit.penalty= conquer::conquer.cv.reg(X=MXZ_ID_SIS, Y=Y, tau = tau_temp, penalty = tolower(penalty))
    beta.penalty = fit.penalty$coeff.min[2:(d+1)]
    
    #---------- Step 3: Mediator significance testing 
    message("        Step 3: Joint significance test ...", "     (", format(Sys.time(), "%X"), ")")
    
    beta_fit_penalty <- matrix(0,1,p)
    ID_penalty <- ID_SIS[which(beta.penalty != 0)]     # the index of nonzero mediators
    MXZ_penalty <- cbind(M[,ID_penalty],XZ)            # cbind M[ID_penalty], X, and Z
    
    fit_rq_penalty <- quantreg::rq(Y ~ MXZ_penalty, tau = tau_temp, method="fn", model=TRUE)
    rq_est_penalty <- summary(fit_rq_penalty, covariance = TRUE, se = "boot")$coefficients
    beta_hat_penalty <- matrix(rq_est_penalty[2:(length(ID_penalty)+1),1])
    beta_SE_penalty  <- matrix(rq_est_penalty[2:(length(ID_penalty)+1),2])
    beta_fit_penalty[ID_penalty] <-  beta_hat_penalty
    
    #--- the p-values
    P_raw_2k <- 2*(1-pnorm(abs(alpha_est)/alpha_SE,0,1)) # the p-values of \alpha
    P_raw_1k_penalty  <-  2*(1-pnorm(abs(beta_hat_penalty)/beta_SE_penalty,0,1))
    P_max_k_penalty <- apply(cbind(P_raw_1k_penalty, t(t(P_raw_2k[ID_penalty]))),1,max)
    
    #-- the Pmax method
    P_penalty_Pmax  <- P_max_k_penalty*length(ID_penalty)
    sig_ind <- which(P_penalty_Pmax < Bonfcut)
    ID_Non_penalty_Pmax <- ID_penalty[sig_ind] # the ID of significant M by JS
    
    if(length(ID_Non_penalty_Pmax) > 0)
    {
      out_result <- rbind(out_result, 
                          data.frame(ID = M_ID_name[ID_Non_penalty_Pmax], 
                                     alpha = alpha_est[ID_Non_penalty_Pmax], 
                                     alpha_se = alpha_SE[ID_Non_penalty_Pmax], 
                                     beta = beta_hat_penalty[sig_ind], 
                                     beta_se = beta_SE_penalty[sig_ind],
                                     Bonferroni.p = P_max_k_penalty[sig_ind],
                                     tau = tau_temp))
      if(verbose) message(paste0("        ", length(ID_Non_penalty_Pmax), " significant mediator(s) identified."))
    } else {
      if(verbose) message("        No significant mediator identified.")
    }
  }
  
  message("Done!", "     (", format(Sys.time(), "%X"), ")")
  return(out_result)
}
