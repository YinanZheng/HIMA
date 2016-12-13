#' High-dimensional Mediation Analysis
#' 
#' \code{hima} is used to estimate and test high-dimensional mediation effects.
#' 
#' @param X a vector of exposure. 
#' @param Y a vector of outcome. Can be either continuous or binary (0-1).
#' @param M a data frame or matrix of high-dimensional mediators. Rows represent samples, columns 
#' represent variables.
#' @param COV a data frame or matrix of covariates dataset. Default = \code{NULL}.
#' @param family either 'gaussian' or 'binomial', depending on the data type of outcome (\code{Y}). See 
#' \code{\link{ncvreg}}
#' @param penalty the penalty to be applied to the model. Either 'MCP' (the default), 'SCAD', or 
#' 'lasso'. See \code{\link{ncvreg}}.
#' @param topN an integer specifying the number of top markers from sure independent screening. 
#' Default = \code{NULL}. If \code{NULL}, \code{ceiling(n/log10(n))} will be used as \code{topN}, 
#' where \code{n} is the sample size. 
#' @param parallel logical. Enable parallel computing feature? Default = \code{TRUE}.
#' @param ncore number of cores to run parallel computation. Valid when \code{parallel == TRUE}. 
#' By default max number of cores available in the machine will be utilized. If \code{ncore = 1}, 
#' no parallel computating is allowed.
#' @param verbose logical. Should the function be verbose?
#' @param ... other arguments passed to \code{\link{ncvreg}}.
#' 
#' @return A data.frame containing mediation testing results of selected mediators.
#'
#' @examples
#' n <- 300  # sample size
#' p <- 10000 # the dimension of covariates
#' 
#' # the regression coefficients alpha (exposure --> mediators)
#' alpha <- rep(0, p) 
#' 
#' # the regression coefficients beta (mediators --> outcome)
#' beta1 <- rep(0, p) # for continuous outcome
#' beta2 <- rep(0, p) # for binary outcome

#' # the first four markers are true mediators
#' alpha[1:4] <- c(0.45,0.5,0.55,0.6)
#' beta1[1:4] <- c(0.5,0.45,0.4,0.35)
#' beta2[1:4] <- c(1.5,1.45,1.4,1.35)
#'
#' # these are also not true mediators
#' alpha[7:8] <- 0.5
#' beta1[5:6] <- 0.4
#' beta2[5:6] <- 1.4
#' 
#' # Generate simulation data
#' simdat_cont = simHIMA(n, p, alpha, beta1, seed=2016) 
#' simdat_bin = simHIMA(n, p, alpha, beta2, binaryOutcome = TRUE, seed=2016) 

#' # Run HIMA with MCP penalty by default
#' hima.fit <- hima(simdat_cont$X, simdat_cont$Y, simdat_cont$M, 
#' parallel = FALSE) 
#' head(hima.fit)
#' 
#' hima.logistic.fit <- hima(simdat_bin$X, simdat_bin$Y, simdat_bin$M, 
#' family = "binomial", , parallel = FALSE) 
#' head(hima.logistic.fit)
#' 
#' @export
hima <- function(X, Y, M, COV = NULL, 
                 family = c("gaussian", "binomial"), 
                 penalty = c("MCP", "SCAD", "lasso"), 
                 topN = NULL, 
                 parallel = TRUE, 
                 ncore = NULL, 
                 verbose = TRUE, 
                 ...) {
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  
  if (is.null(ncore)) 
    ncore <- parallel::detectCores()
  
  n <- nrow(M)
  p <- ncol(M)
  if (is.null(topN)) {
    if (family == "binomial") d <- ceiling(n/(2*log(n))) else d <- ceiling(2 * n/log(n)) 
  } else {
    d <- topN  # the number of top mediators that associated with exposure (X)
  }
  
  if(verbose) message("Step 1: Screening...", "     (", Sys.time(), ")")
  
  if(family == "binomial")
  {
    alpha = SIS_alpha <- himasis(NA, M, X, COV, glm.family = "gaussian", modelstatement = "Mone ~ X", 
                             parallel = parallel, ncore = ncore, verbose)
    SIS_p <- SIS_alpha[2,]
  } else {
    SIS_beta <- himasis(Y, M, X, COV, glm.family = family, modelstatement = "Y ~ Mone + X", 
                             parallel = parallel, ncore = ncore, verbose)
    SIS_p <- SIS_beta[2,]
  }

  SIS_p_sort <- sort(SIS_p)
  ID <- which(SIS_p <= SIS_p_sort[d])  # the index of top mediators
  M_SIS <- M[, ID]
  XM <- cbind(M_SIS, X)
  
  if(verbose) message("Step 2: Penalized estimate (", penalty, ") ...", "     (", 
          Sys.time(), ")")
  
  if (is.null(COV)) {
    fit <- ncvreg(XM, Y, family = family, 
                  penalty = penalty, 
                  penalty.factor = c(rep(1, ncol(M_SIS)), 0))
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1]
    conf.names <- colnames(COV)
    XM_COV <- cbind(XM, COV)
    fit <- ncvreg(XM_COV, Y, family = family, 
                  penalty = penalty, 
                  penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(COV))), ...)
  }
  # plot(fit)
  
  lam <- fit$lambda[which.min(BIC(fit))]
  Coefficients <- coef(fit, lambda = lam)
  est <- Coefficients[2:(d + 1)]
  ID_1_non <- which(est != 0)
  beta_est <- est[ID_1_non]  # The non-zero MCP estimators of \beta
  ID_test <- ID[ID_1_non]  # The index of the ID of non-zero beta in the regression of Y= c + rX + b_1 M_1 + ...+ b_p M_p
  ## 
  
  if(verbose) message("Step 3: Joint significance test ...", 
                      "     (", Sys.time(), ")")
  
  if(family == "binomial")
  {
    alpha <- alpha[,ID_test]
    # modelstatement <- eval(parse(text = paste0("Y ~ ", paste0(names(ID_test), collapse = "+"))))
    # if (is.null(COV)) {
    #   datarun <- data.frame(Y = Y, M[,ID_test])
    #   modelstatement <- modelstatement
    # } else {
    #   datarun <- data.frame(Y = Y, M[,ID_test], COV = COV)
    #   modelstatement <- eval(parse(text = (paste0(modelstatement, "+", 
    #                                               paste0(paste0("COV.", conf.names), collapse = "+")))))
    # }
    # 
    # beta <- glm(modelstatement, data = datarun, family = family)
    # 
  } else {
    alpha <- himasis(NA, M[, ID_test], X, COV, glm.family = "gaussian", 
                     modelstatement = "Mone ~ X", parallel = parallel, ncore = ncore, 
                     verbose = FALSE)
  }

  
  alpha_est_ID_test <- as.numeric(alpha[1, ])  #  the estimator for alpha
  P_adjust_alpha <- length(ID_test) * alpha[2, ]  # The adjusted p-value for alpha
  P_adjust_alpha[P_adjust_alpha > 1] <- 1
  
  alpha_est <- alpha_est_ID_test
  
  #################### the following is the post-test based on the oracle property of the
  #################### MCP penalty
  if (is.null(COV)) {
    YMX <- data.frame(Y = Y, M[, ID_test], X = X)
  } else {
    YMX <- data.frame(Y = Y, M[, ID_test], X = X, COV)
  }
  
  res <- summary(glm(Y ~ ., family = family, data = YMX))$coefficients
  est <- res[2:(length(ID_test) + 1), 1]  # The estimator for alpha
  P_adjust_beta <- length(ID_test) * res[2:(length(ID_test) + 1), 4]  # The adjused p-value for beta
  P_adjust_beta[P_adjust_beta > 1] <- 1
  
  ab_est <- alpha_est * beta_est
  ## 
  PA <- rbind(P_adjust_beta, P_adjust_alpha)
  P_value <- apply(PA, 2, max)
  
  P_value_0 <- P_value_1 <- rep(1, p)
  P_value_0[ID_test] <- P_value
  
  PB <- rbind(P_value_0, P_value_1)
  P_value <- apply(PB, 2, min)
  
  
  # Total effect
  if (is.null(COV)) {
    YX <- data.frame(Y = Y, X = X)
  } else {
    YX <- data.frame(Y = Y, X = X, COV)
  }
  
  gamma_est <- coef(glm(Y ~ ., family = family, data = YX))[2]
  
  results <- data.frame(alpha = alpha_est, beta = beta_est, gamma = gamma_est, 
                        `alpha*beta` = ab_est, `% total effect` = ab_est/gamma_est * 100, 
                        `p-value` = P_value[ID_test], check.names = FALSE)
  
  if(verbose) message("Done!", "     (", Sys.time(), ")")
  
  closeAllConnections()
  registerDoSEQ()
  
  return(results)
}