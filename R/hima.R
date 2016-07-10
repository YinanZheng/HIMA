#' High-dimensional Mediation Analysis
#' 
#' \code{hima} is used to estimate and test high-dimensional mediation effects.
#' 
#' @param X a vector of exposure. 
#' @param Y a vector of outcome.
#' @param M a data frame or matrix of high-dimensional mediators. Rows represent samples, columns represent variables.
#' @param COV a data frame or matrix of covariates dataset. Default = \code{NULL}.
#' @param family Either "gaussian", "binomial", or "poisson", depending on the response. See \code{\link{ncvreg}}
#' @param penalty	The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso". See \code{\link{ncvreg}}.
#' @param topN an integer specifying the number of top markers from sure independent screening. Default = \code{NULL}. If \code{NULL}, \code{ceiling(n/log10(n))} will be used as \code{topN}, where \code{n} is the sample size. 
#' @param scale logical. Whether standardization should be applied to the high dimensional mediator data \code{M}. Default = \code{TRUE}.
#' @param parallel logical. Enable parallel computing feature? Default = \code{TRUE}.
#' @param ncore number of cores to run parallel. Effective when \code{parallel = TRUE}. Default = \code{2}.
#' 
#' @return a data.frame contains mediation testing results of selected mediators.
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
#' simdat = simHIMA(n, p, alpha, beta, seed=2016) # Generate simulation data
#' 
#' hima(simdat$X, simdat$Y, simdat$M, penalty = "SCAD", ncore=detectCores()) # Run HIMA with MCP penalty
#' hima(simdat$X, simdat$Y, simdat$M, penalty = "SCAD", ncore=detectCores()) # Run HIMA with SCAD penalty
#' hima(simdat$X, simdat$Y, simdat$M, penalty = "lasso", ncore=detectCores()) # Run HIMA with lasso penalty
#' 
#' @export
hima <- function(X, Y, M, COV = NULL, family="gaussian", penalty="MCP", topN=NULL, scale = TRUE, parallel=TRUE, ncore=2)
{
  n = nrow(M)
  p = ncol(M)
  if(is.null(topN)) d <- ceiling(2*n/log10(n)) else d = topN  # the number of top mediators that associated with exposure (X)
  
  cat(paste0("Step 1: Screening...", "     (",Sys.time(),")"),"\n")
  SIS_P_value_b = himasis(Y, M, X, COV, modelstatement = "Y ~ Mone + X", scale = scale, parallel = TRUE, ncore = ncore)[2,]
  SIS_P_value_1 = sort(SIS_P_value_b)
  ID <- which(SIS_P_value_b <= SIS_P_value_1[d]) # the index of top mediators
  M_SIS <- M[,ID]
  XM <- cbind(M_SIS, X)
  
  cat(paste0("Step 2: Penalized estimate (", penalty,") ...", "     (",Sys.time(),")"),"\n")
  
  if(is.null(COV)){
    fit <- ncvreg(XM, Y, family = family, penalty = penalty, penalty.factor=c(rep(1, ncol(M_SIS)), 0))
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~.,COV))[,-1]
    XM_COV = cbind(XM, COV)
    fit <- ncvreg(XM_COV, Y, family = family, penalty = penalty, penalty.factor=c(rep(1, ncol(M_SIS)), rep(0, 1+ncol(COV))))
  }

  lam <- fit$lambda[which.min(BIC(fit))]
  Coefficients <- coef(fit, lambda=lam)
  est <- Coefficients[2:(d+1)]
  ID_1_non <- which(est !=0)
  beta_est  <- est[ID_1_non]    # The non-zero MCP estimators of \beta
  ID_test <- ID[ID_1_non]               # The index of the ID of non-zero beta in the regression of Y= c + rX + b_1 M_1 + ...+ b_p M_p
  ##

  cat(paste0("Step 3: Joint significance test ...", "     (",Sys.time(),")"),"\n")
  
  alpha = himasis(NA, M[,ID_test], X, COV, modelstatement = "Mone ~ X", scale = scale, parallel = TRUE, ncore = ncore)
  
  alpha_est_ID_test <- as.numeric(alpha[1,]) #  the estimator for alpha
  P_adjust_alpha <- length(ID_test) * alpha[2,] # The adjusted p-value for alpha
  P_adjust_alpha[P_adjust_alpha>1] = 1
  
  alpha_est  <- alpha_est_ID_test
  
  ####################   the following is the post-test based on the oracle property of the MCP penalty
  if(is.null(COV)){
    YMX <- data.frame(Y = Y, M[,ID_test], X = X)
  } else {
    YMX <- data.frame(Y = Y, M[,ID_test], X = X, COV)
  }
  
  res = summary(lm(Y~., data = YMX))$coefficients
  est <-  res[2:(length(ID_test)+1),1] # The estimator for alpha
  P_adjust_beta <- length(ID_test) * res[2:(length(ID_test)+1),4] # The adjused p-value for beta
  P_adjust_beta[P_adjust_beta>1] = 1
  
  ab_est <- alpha_est * beta_est
  ##
  PA <- rbind(P_adjust_beta,P_adjust_alpha)
  P_value <- apply(PA,2,max)   
  
  P_value_0 = P_value_1 <- rep(1, p)
  P_value_0[ID_test] <- P_value
  
  PB <- rbind(P_value_0,P_value_1)
  P_value <- apply(PB,2,min)
  

  # Total effect
  if(is.null(COV)){
    YX <- data.frame(Y = Y, X = X)
  } else {
    YX <- data.frame(Y = Y, X = X, COV)
  }
  
  gamma_est = coef(lm(Y~., data = YX))[2]

  results = data.frame(alpha = alpha_est, 
                       beta = beta_est, 
                       gamma = gamma_est, 
                       'alpha*beta' = ab_est, 
                       '% total effect' = ab_est/gamma_est*100,
                       'p-value' = P_value[ID_test], check.names = FALSE)
  
  cat(paste0("Done!", "     (",Sys.time(),")"),"\n")
  return(results)
}