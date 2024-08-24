#' This is the function for high-dimensional mediation analysis using de-biased lasso

#' HIMA with de-biased lasso
#' 
#' \code{dblassoHIMA} is used to estimate and test high-dimensional mediation effects using de-biased lasso penalty.
#' 
#' @param X a vector of exposure. 
#' @param Y a vector of outcome. Can be either continuous or binary (0-1).
#' @param M a \code{data.frame} or \code{matrix} of high-dimensional mediators. Rows represent samples, columns represent variables.
#' @param COV a \code{data.frame} or \code{matrix} of covariates dataset for testing the association M ~ X and Y ~ M.
#' @param Y.family either 'gaussian' (default) or 'binomial', depending on the data type of outcome (\code{Y}). This parameter is passed 
#' to function \code{lasso.proj} in R package \code{hdi} for de-biased lasso penalization.
#' @param topN an integer specifying the number of top markers from sure independent screening. 
#' Default = \code{NULL}. If \code{NULL}, \code{topN} will be either \code{ceiling(n/log(n))} if 
#' \code{Y.family = 'gaussian'}, or \code{ceiling(n/(2*log(n)))} if \code{Y.family = 'binomial'}, 
#' where \code{n} is the sample size. If the sample size is greater than topN (pre-specified or calculated), all 
#' mediators will be included in the test (i.e. low-dimensional scenario).
#' @param scale logical. Should the function scale the data? Default = \code{TRUE}.
#' @param FDRcut HDMT pointwise FDR cutoff applied to select significant mediators. Default = \code{0.05}. 
#' @param verbose logical. Should the function be verbose? Default = \code{FALSE}.
#' 
#' @return A data.frame containing mediation testing results of significant mediators (FDR <\code{FDRcut}). 
#' \describe{
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
#' @references Perera C, Zhang H, Zheng Y, Hou L, Qu A, Zheng C, Xie K, Liu L. 
#' HIMA2: high-dimensional mediation analysis and its application in epigenome-wide DNA methylation data. 
#' BMC Bioinformatics. 2022. DOI: 10.1186/s12859-022-04748-1. PMID: 35879655; PMCID: PMC9310002
#' 
#' @examples
#' \dontrun{
#' # Note: In the following examples, M1, M2, and M3 are true mediators.
#' data(himaDat)
#' 
#' # When Y is continuous and normally distributed
#' # Example 1 (continuous outcome): 
#' head(himaDat$Example1$PhenoData)
#' 
#' dblassohima.fit <- dblassoHIMA(X = himaDat$Example1$PhenoData$Treatment, 
#'                                Y = himaDat$Example1$PhenoData$Outcome, 
#'                                M = himaDat$Example1$Mediator,
#'                                COV = himaDat$Example1$PhenoData[, c("Sex", "Age")],
#'                                Y.family = 'gaussian',
#'                                scale = FALSE, # Disabled only for simulation data
#'                                verbose = TRUE) 
#' dblassohima.fit
#' 
#' # When Y is binary (should specify Y.family)
#' # Example 2 (binary outcome): 
#' head(himaDat$Example2$PhenoData)
#' 
#' dblassohima.logistic.fit <- dblassoHIMA(X = himaDat$Example2$PhenoData$Treatment,
#'                                         Y = himaDat$Example2$PhenoData$Disease,
#'                                         M = himaDat$Example2$Mediator,
#'                                         COV = himaDat$Example2$PhenoData[, c("Sex", "Age")],
#'                                         Y.family = 'binomial',
#'                                         scale = FALSE, # Disabled only for simulation data
#'                                         verbose = TRUE)
#' dblassohima.logistic.fit
#' }
#' 
#' @export
dblassoHIMA<-function(X, Y, M, COV = NULL, 
                      Y.family = c("gaussian", "binomial"),
                      topN = NULL, 
                      scale = TRUE,
                      FDRcut = 0.05,
                      verbose = FALSE)
{
  Y.family <- match.arg(Y.family)
  
  n <- nrow(M)
  p <- ncol(M)
  
  if(scale)
  {
    X <- scale(X)
    M <- scale(M)
    if(!is.null(COV)) COV <- scale(COV)
    if(verbose) message("Data scaling is completed.")
  } else {
    X <- as.matrix(X)
    M <- as.matrix(M)
    if(!is.null(COV)) COV <- as.matrix(COV)
  }
  
  if(is.null(COV)) {MZX <- cbind(M,X); XZ <- X; q <- 0} else {MZX <- cbind(M,COV,X); XZ <- cbind(X, COV); q <- ncol(COV)}
  
  if(is.null(topN)) {
    d <- ceiling(2 * n/log(n)) 
  } else {
    d <- topN  # the number of top mediators that associated with exposure (X)
  }
  
  d <- min(p, d) # if d > p select all mediators
  
  M_ID_name <- colnames(M)
  if(is.null(M_ID_name)) M_ID_name <- seq_len(p)
  
  #########################################################################
  ########################### (Step 1) SIS step ########################### 
  #########################################################################
  message("Step 1: Sure Independent Screening ...", "  (", format(Sys.time(), "%X"), ")")
  
  # the number of top mediators that associated with exposure (X)
  if(is.null(topN)) d_0 <- ceiling(2 * n/log(n)) else d_0 <- topN  
  
  d_0 <- min(p, d_0) # if d > p select all mediators
  
  beta_SIS <- matrix(0,1,p) 
  
  # Estimate the regression coefficients beta (mediators --> outcome)
  for (i in 1:p){
    ID_S <- c(i, (p+1):(p+q+1))
    MZX_SIS <- MZX[,ID_S]
    fit <- lsfit(MZX_SIS,Y,intercept = TRUE)
    beta_SIS[i] <- fit$coefficients[2]
  }
  
  # Estimate the regression coefficients alpha (exposure --> mediators)
  alpha_SIS <- matrix(0,1,p)
  for (i in 1:p){
    fit_a  <- lsfit(XZ,M[,i],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    alpha_SIS[i] <- est_a
  }
  
  # Select the d_0 number of mediators with top largest effect 
  ab_SIS <- alpha_SIS*beta_SIS
  ID_SIS  <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[d_0])
  d <- length(ID_SIS)
  
  if(verbose) message("        Top ", d, " mediators are selected: ", paste0(M_ID_name[ID_SIS], collapse = ", "))
  
  #########################################################################
  ################### (Step 2) De-biased Lasso Estimates ##################
  #########################################################################
  message("Step 2: De-biased Lasso Estimates ...", "   (", format(Sys.time(), "%X"), ")")
  
  if(verbose)
  {
    if(is.null(COV)) 
    {message("        No covariate was adjusted.")} 
    else
    {message("        Adjusting for covariate(s): ", paste0(colnames(COV), collapse = ", "))}
  }
  
  P_beta_SIS <- matrix(0,1,d)
  beta_DLASSO_SIS_est <- matrix(0,1,d)
  beta_DLASSO_SIS_SE <- matrix(0,1,d)
  MZX_SIS <- MZX[,c(ID_SIS, (p+1):(p+q+1))]
  
  DLASSO_fit <- suppressMessages(hdi::lasso.proj(x=MZX_SIS, y=Y, family = Y.family))
  beta_DLASSO_SIS_est <- DLASSO_fit$bhat[1:d]
  beta_DLASSO_SIS_SE <- DLASSO_fit$se
  P_beta_SIS <- t(DLASSO_fit$pval[1:d])
  
  ################### Estimate alpha ################
  alpha_SIS_est <- matrix(0,1,d)
  alpha_SIS_SE <- matrix(0,1,d)
  P_alpha_SIS <- matrix(0,1,d)
  
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
  ################ (step 3) The multiple-testing  procedure ###############
  #########################################################################
  message("Step 3: Joint significance test ...", "     (", format(Sys.time(), "%X"), ")")
  
  PA <- cbind(t(P_alpha_SIS),(t(P_beta_SIS)))
  P_value <- apply(PA,1,max)  #The joint p-values for SIS variable
  
  N0 <- dim(PA)[1]*dim(PA)[2]
  input_pvalues <- PA + matrix(runif(N0,0,10^{-10}),dim(PA)[1],2)
  
  # Estimate the proportions of the three component nulls
  nullprop <- null_estimation(input_pvalues)
  
  # Compute the estimated pointwise FDR for every observed p-max
  fdrcut  <- HDMT::fdr_est(nullprop$alpha00,
                           nullprop$alpha01,
                           nullprop$alpha10, 
                           nullprop$alpha1,
                           nullprop$alpha2,
                           input_pvalues,exact=0)
  
  ID_fdr <- which(fdrcut <= FDRcut)
  
  # Following codes extract the estimates for mediators with fdrcut<=0.05
  beta_hat_est <- beta_DLASSO_SIS_est[ID_fdr]
  beta_hat_SE  <- beta_DLASSO_SIS_SE[ID_fdr]
  
  alpha_hat_est <-  alpha_SIS_est[ID_fdr]
  alpha_hat_SE  <-  alpha_SIS_SE[ID_fdr]
  
  P.value_raw <- P_value[ID_fdr]
  
  # Indirect effect
  IDE <- beta_hat_est*alpha_hat_est # mediation(indirect) effect
  
  # # Total effect
  # if(is.null(COV)) {
  #   YX <- data.frame(Y = Y, X = X)
  # } else {
  #   YX <- data.frame(Y = Y, X = X, COV)
  # }
  # 
  # gamma_est <- coef(glm(Y ~ ., family = Y.family, data = YX))[2]
  
  if(length(ID_fdr) > 0)
  {
    results <- data.frame(Index = M_ID_name[ID_fdr], 
                          alpha_hat = alpha_hat_est,
                          alpha_se = alpha_hat_SE,
                          beta_hat = beta_hat_est, 
                          beta_se = beta_hat_SE,
                          IDE = IDE, 
                          rimp = abs(IDE)/sum(abs(IDE)) * 100, 
                          pmax = P.value_raw, check.names = FALSE)
    
    if(verbose) message(paste0("        ", length(ID_fdr), " significant mediator(s) identified."))
  } else {
    if(verbose) message("        No significant mediator identified.")
    results = NULL
  }
  
  message("Done!", "     (", format(Sys.time(), "%X"), ")")
  
  return(results)
}