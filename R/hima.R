#' High-dimensional Mediation Analysis
#' 
#' \code{hima} is used to estimate and test high-dimensional mediation effects.
#' 
#' @param X a vector of exposure. 
#' @param Y a vector of outcome. Can be either continuous or binary (0-1).
#' @param M a \code{data.frame} or \code{matrix} of high-dimensional mediators. Rows represent samples, columns 
#' represent variables.
#' @param COV.XM a \code{data.frame} or \code{matrix} of covariates dataset for testing the association \code{M ~ X}. 
#' Covariates specified here will not participate penalization. Default = \code{NULL}. If the covariates 
#' contain mixed data types, please make sure all categorical variables are properly formatted as \code{factor} type.
#' @param COV.MY a \code{data.frame} or \code{matrix} of covariates dataset for testing the association \code{Y ~ M}. 
#' Covariates specified here will not participate penalization. If not specified, the same set of covariates for 
#' \code{M ~ X} will be applied. Using different sets of covariates is allowed but this needs to be handled carefully.
#' @param family either 'gaussian' or 'binomial', depending on the data type of outcome (\code{Y}). See 
#' \code{\link{ncvreg}}
#' @param penalty the penalty to be applied to the model. Either 'MCP' (the default), 'SCAD', or 
#' 'lasso'. See \code{\link{ncvreg}}.
#' @param topN an integer specifying the number of top markers from sure independent screening. 
#' Default = \code{NULL}. If \code{NULL}, \code{topN} will be either \code{ceiling(n/log(n))} if 
#' \code{family = 'gaussian'}, or \code{ceiling(n/(2*log(n)))} if \code{family = 'binomial'}, 
#' where \code{n} is the sample size. If the sample size is greater than topN (pre-specified or calculated), all 
#' mediators will be included in the test (i.e. low-dimensional scenario).
#' @param parallel logical. Enable parallel computing feature? Default = \code{TRUE}.
#' @param ncore number of cores to run parallel computing Valid when \code{parallel == TRUE}. 
#' By default max number of cores available in the machine will be utilized.
#' @param verbose logical. Should the function be verbose? Default = \code{FALSE}.
#' @param ... other arguments passed to \code{\link{ncvreg}}.
#' 
#' @return A data.frame containing mediation testing results of selected mediators. 
#' \itemize{
#'     \item{alpha: }{coefficient estimates of exposure (X) --> mediators (M).}
#'     \item{beta: }{coefficient estimates of mediators (M) --> outcome (Y) (adjusted for exposure).}
#'     \item{gamma: }{coefficient estimates of exposure (X) --> outcome (Y) (total effect).}
#'     \item{alpha*beta: }{mediation effect.}
#'     \item{\% total effect: }{alpha*beta / gamma. Percentage of the mediation effect out of the total effect.}
#'     \item{Bonferroni.p: }{statistical significance of the mediator (Bonferroni procedure).}
#'     \item{BH.FDR: }{statistical significance of the mediator (Benjamini-Hochberg procedure).}
#' }
#' 
#' @references Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J, Just A, Colicino E, Vokonas P, Zhao L, 
#' Lv J, Baccarelli A, Hou L, Liu L. Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies. 
#' Bioinformatics. 2016. DOI: 10.1093/bioinformatics/btw351. PMID: 27357171. PMCID: PMC5048064
#' 
#' @examples
#' n <- 200  # sample size
#' p <- 200 # the dimension of covariates
#' 
#' # the regression coefficients alpha (exposure --> mediators)
#' alpha <- rep(0, p) 
#' 
#' # the regression coefficients beta (mediators --> outcome)
#' beta1 <- rep(0, p) # for continuous outcome
#' beta2 <- rep(0, p) # for binary outcome
#' 
#' # the first four markers are true mediators
#' alpha[1:4] <- c(0.45,0.5,0.6,0.7)
#' beta1[1:4] <- c(0.55,0.6,0.65,0.7)
#' beta2[1:4] <- c(1.45,1.5,1.55,1.6)
#'
#' # these are not true mediators
#' alpha[7:8] <- 0.5
#' beta1[5:6] <- 0.8
#' beta2[5:6] <- 1.7
#' 
#' # Generate simulation data
#' simdat_cont = simHIMA(n, p, alpha, beta1, seed=1029) 
#' simdat_bin = simHIMA(n, p, alpha, beta2, binaryOutcome = TRUE, seed=1029) 
#' 
#' # Run HIMA with MCP penalty by default
#' # When Y is continuous (default)
#' hima.fit <- hima(simdat_cont$X, simdat_cont$Y, simdat_cont$M, verbose = TRUE) 
#' hima.fit
#' 
#' # When Y is binary (should specify family)
#' hima.logistic.fit <- hima(simdat_bin$X, simdat_bin$Y, simdat_bin$M, 
#' family = "binomial", verbose = TRUE) 
#' hima.logistic.fit
#' 
#' @export
hima <- function(X, Y, M, COV.XM = NULL, COV.MY = COV.XM, 
                 family = c("gaussian", "binomial"), 
                 penalty = c("MCP", "SCAD", "lasso"), 
                 topN = NULL, 
                 parallel = FALSE, 
                 ncore = 1, 
                 verbose = FALSE, 
                 ...) {
    family <- match.arg(family)
    penalty <- match.arg(penalty)
    
    if (parallel & (ncore == 1)) ncore <- parallel::detectCores()

    n <- nrow(M)
    p <- ncol(M)
    
    if (is.null(topN)) {
      if (family == "binomial") d <- ceiling(n/(2*log(n))) else d <- ceiling(2 * n/log(n)) 
    } else {
      d <- topN  # the number of top mediators that associated with exposure (X)
    }
  
    d <- min(p, d) # if d > p select all mediators
    
    #########################################################################
    ################################ STEP 1 #################################
    #########################################################################
    message("Step 1: Sure Independent Screening ...", "     (", Sys.time(), ")")
    
    if(family == "binomial")
    {
      # Screen M using X given the limited information provided by Y (binary)
      # Therefore the family is still gaussian
      if(verbose) message("    Screening M using the association between X and M: ", appendLF = FALSE)
      alpha = SIS_Results <- himasis(NA, M, X, COV.XM, glm.family = "gaussian", modelstatement = "Mone ~ X", 
                               parallel = parallel, ncore = ncore, verbose, tag = "Sure Independent Screening")
      SIS_Pvalue <- SIS_Results[2,]
    } else if (family == "gaussian"){
      # Screen M using Y (continuous)
      if(verbose) message("    Screening M using the association between M and Y: ", appendLF = FALSE)
      SIS_Results <- himasis(Y, M, X, COV.MY, glm.family = family, modelstatement = "Y ~ Mone + X", 
                               parallel = parallel, ncore = ncore, verbose, tag = "Sure Independent Screening")
      SIS_Pvalue <- SIS_Results[2,]
    } else {
      stop(paste0("Family ", family, " is not supported."))
    }
    # Note: ranking using p on un-standardized data is equivalent to ranking using beta on standardized data
    SIS_Pvalue_sort <- sort(SIS_Pvalue)
    ID <- which(SIS_Pvalue <= SIS_Pvalue_sort[d])  # the index of top mediators
    if(verbose) message("    Top ", length(ID), " mediators are selected: ", paste0(names(SIS_Pvalue_sort[seq_len(d)]), collapse = ","))
  
    M_SIS <- M[, ID]
    XM <- cbind(M_SIS, X)
    
    #########################################################################
    ################################ STEP 2 #################################
    #########################################################################
    message("Step 2: Penalized estimate (", penalty, ") ...", "     (", Sys.time(), ")")
    
    ## Based on the screening results in step 1. We will find the most influential M on Y.
    if (is.null(COV.MY)) {
      fit <- ncvreg(XM, Y, family = family, 
                    penalty = penalty, 
                    penalty.factor = c(rep(1, ncol(M_SIS)), 0), ...)
    } else {
      COV.MY <- data.frame(COV.MY)
      COV.MY <- data.frame(model.matrix(~., COV.MY))[, -1]
      conf.names <- colnames(COV.MY)
      if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
      XM_COV <- cbind(XM, COV.MY)
      fit <- ncvreg(XM_COV, Y, family = family, 
                    penalty = penalty, 
                    penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(COV.MY))), ...)
    }
    # plot(fit)
    
    lam <- fit$lambda[which.min(BIC(fit))]
    if(verbose) message("    Tuning parameter lambda selected: ", lam)
    Coefficients <- coef(fit, lambda = lam)
    est <- Coefficients[2:(d + 1)]
    ID_1_non <- which(est != 0)
    if(length(ID_1_non) == 0)
    {
      if(verbose) message("    All ", penalty, " beta estimates of the ", length(ID), " mediators are zero. Please check ncvreg package for more options (?ncvreg).")
    } else {
    if(verbose) message("    Non-zero ", penalty, " beta estimate(s) of mediator(s) found: ", paste0(names(ID_1_non), collapse = ","))
    beta_est <- est[ID_1_non]  # The non-zero MCP estimators of beta
    ID_test <- ID[ID_1_non]  # The index of the ID of non-zero beta in Y ~ M
    ## 
    
    if(family == "binomial")
    {
      ## This has been done in step 1 (when Y is binary)
      alpha <- alpha[,ID_test, drop = FALSE]
    } else {
      if(verbose) message("    Estimating alpha (effect of X on M): ", appendLF = FALSE)
      alpha <- himasis(NA, M[, ID_test, drop = FALSE], X, COV.XM, glm.family = "gaussian", 
                       modelstatement = "Mone ~ X", parallel = FALSE, ncore = ncore, 
                       verbose, tag = "site-by-site ordinary least squares estimation")
    }
    
    #########################################################################
    ################################ STEP 3 #################################
    #########################################################################
    message("Step 3: Joint significance test ...", "     (", Sys.time(), ")")
    
    alpha_est_ID_test <- as.numeric(alpha[1, ])  #  the estimator for alpha
    P_adjust_alpha <- length(ID_test) * alpha[2, ]  # the adjusted p-value for alpha (bonferroni)
    P_adjust_alpha[P_adjust_alpha > 1] <- 1
    P_fdr_alpha <- p.adjust(alpha[2, ], "fdr")  # the adjusted p-value for alpha (FDR)
    
    alpha_est <- alpha_est_ID_test
    
    ## Post-test based on the oracle property of the MCP penalty
    if (is.null(COV.MY)) {
      YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X)
    } else {
      YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X, COV.MY)
    }
    
    res <- summary(glm(Y ~ ., family = family, data = YMX))$coefficients
    est <- res[2:(length(ID_test) + 1), 1]  # the estimator for beta
    P_adjust_beta <- length(ID_test) * res[2:(length(ID_test) + 1), 4]  # the adjused p-value for beta (bonferroni)
    P_adjust_beta[P_adjust_beta > 1] <- 1
    P_fdr_beta <- p.adjust(res[2:(length(ID_test) + 1), 4], "fdr")  # the adjusted p-value for beta (FDR)
    
    ab_est <- alpha_est * beta_est
    
    ## Use the maximum value as p value 
    PA <- rbind(P_adjust_beta, P_adjust_alpha)
    P_value <- apply(PA, 2, max)
    
    FDRA <- rbind(P_fdr_beta, P_fdr_alpha)
    FDR <- apply(FDRA, 2, max)
    
    # Total effect
    if (is.null(COV.MY)) {
      YX <- data.frame(Y = Y, X = X)
    } else {
      YX <- data.frame(Y = Y, X = X, COV.MY)
    }
    
    gamma_est <- coef(glm(Y ~ ., family = family, data = YX))[2]
    
    results <- data.frame(alpha = alpha_est, beta = beta_est, gamma = gamma_est, 
                          `alpha*beta` = ab_est, `% total effect` = ab_est/gamma_est * 100, 
                          `Bonferroni.p` = P_value, `BH.FDR` = FDR, check.names = FALSE)
    
    message("Done!", "     (", Sys.time(), ")")
    
    doParallel::stopImplicitCluster()
    
    return(results)
  }
}
