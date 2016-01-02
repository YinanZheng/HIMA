#' High-dimensional Mediation Analysis
#' 
#' \code{hima} is used to estimate and test high-dimensional mediation effects.
#' 
#' @param X a vector of exposure. 
#' @param Y a vector of outcome.
#' @param M a data frame or matrix of high-dimensional mediators. Rows represent samples, columns represent variables.
#' @param COV a data frame or matrix of covariates dataset. Default = \code{NULL}.
#' @param topN an integer specifying the number of top markers from sure independent screening. Default = \code{NULL}. If \code{NULL}, \code{ceiling(n/log10(n))} will be used as \code{topN}, where \code{n} is the sample size. 
#' @param pn number of permutation. Default = \code{500}.
#' @param alpha elastic-net mixing parameter, ranging from \code{0} to \code{1}. alpha = \code{1} is the lasso penalty, and alpha = \code{0} the ridge penalty. Default = \code{0.2}. See \code{\link{glmnet}}.
#' @param parallel logical. Enable parallel computing feature? Default = \code{TRUE}.
#' @param ncore number of cores to run parallel. Effective when \code{parallel = TRUE}. Default = \code{2}.
#' 
#' @return a data.frame contains mediation testing results of mediators that survived from the sure independent screening (see \code{\link{sis}}). Columns include names of mediators (\code{ID}), effects of exposure (\code{a_est_SIS}), effects of mediators (\code{b_est_SIS}), and p-value of mediation test (\code{p}).
#' 
#' @seealso see \code{\link{hima}} using the results from \code{\link{sis}} as input to run HIMA.
#'
#' @examples
#' n <- 100  # sample size
#' p <- 5000 # the dimension of covariates
#' 
#' alpha <- rep(0, p) # the regression coefficients alpha (exposure --> mediators)
#' beta <- rep(0, p) # the regression coefficients beta (mediators --> outcome)
# 
#' alpha[1:3] <- c(0.5, 0.5, 0.3) 
#' beta[1:3] <- c(0.5, 1.2, 0.3) # the first three makers are true mediators.
#' 
#' simdat = simHMA(n, p, alpha, beta, seed=2016) # Generate simulation data
#' 
#' res = hima(simdat$X, simdat$Y, simdat$M, ncore=detectCores()) # Run HIMA
#' head(res)
#' 
#' @export
hima <- function(X, Y, M, COV = NULL, topN=NULL, pn=500, alpha=0.2, parallel=TRUE, ncore=2)
{
  n = nrow(M)
  p = ncol(M)
  if(is.null(topN)) d <- ceiling(n/log10(n)) else d = topN  # the number of top mediators that associated with exposure (X)
  
  ## the following is the SIS idea
  cat(paste0("Start running sure independent screening...", "     (",Sys.time(),")"),"\n")
  a_est = sis(M, X, COV, parallel = parallel, ncore = ncore)
  a_est = a_est[order(abs(a_est), decreasing = T)]
  namd_SIS = names(a_est)[1:d]
  ID <- match(namd_SIS, colnames(M)) # the index of top mediators
  ####
  a_est_SIS <- a_est[1:d]
  M_SIS <- M[,ID]
  X_SIS <- cbind(M_SIS,X)
  
  fit = glmnet(X_SIS, Y, alpha = alpha)
  cv.fit <- cv.glmnet(X_SIS, Y, alpha = alpha)
  Coefficients <- coef(fit, s = cv.fit$lambda.min)
  b_est_SIS <- Coefficients[2:(d+1)]  # the estimator for b
  
  ## the following begin the permutation
  cat(paste0("Start running permutation (repeats=", pn,")...","     (",Sys.time(),")"),"\n")
  PN <- pn
  a_est_per <- matrix(0,PN,d)
  b_est_per <- matrix(0,PN,d)
  for (k in 1: PN){
    if(k %% 100 == 0) cat("-->",k,"times","\n")
    a_est_per[k,] <- sis(M_SIS[sample(1:d,d),], X, COV, parallel = parallel, ncore = ncore)
    fit = glmnet(X_SIS, sample(Y), alpha = alpha)
    cv.fit <- cv.glmnet(X_SIS, sample(Y), alpha = alpha)
    Coefficients <- coef(fit, s = cv.fit$lambda.min)
    b_est_per[k,] <- Coefficients[2:(d+1)]  # the estimator for b
  } # the end of permutation

  P_matrix_a   <- matrix(0,PN,d)
  P_matrix_b <- matrix(0,PN,d)
  for (i in 1:PN){
    P_matrix_a[i,]    <- as.numeric(abs(a_est_per[i,]) >= abs(a_est_SIS))
    P_matrix_b[i,]  <- as.numeric(abs(b_est_per[i,]) >= abs(b_est_SIS))
  }
  
  P_value_a <- colMeans(P_matrix_a)
  P_value_b <- colMeans(P_matrix_b)
  p_mix  <-  rbind(P_value_a,P_value_b)
  P_value <- apply(p_mix,2,max) # the p-value
  
  results = data.frame(ID = namd_SIS, a_est_SIS = a_est_SIS, b_est_SIS = b_est_SIS, p = P_value)
  return(results)
  cat(paste0("Done!", "     (",Sys.time(),")"),"\n")
}