# This is the function classic high-dimensional mediation analysis
#' Classic high-dimensional mediation analysis
#'
#' \code{hima_classic} is used to estimate and test classic high-dimensional mediation effects (linear & logistic regression).
#'
#' @param X a vector of exposure. Do not use \code{data.frame} or \code{matrix}.
#' @param M a \code{data.frame} or \code{matrix} of high-dimensional mediators. Rows represent samples, columns
#' represent variables.
#' @param Y a vector of outcome. Can be either continuous or binary (0-1). Do not use \code{data.frame} or \code{matrix}.
#' @param COV.XM a \code{data.frame} or \code{matrix} of covariates dataset for testing the association \code{M ~ X}.
#' Covariates specified here will not participate penalization. Default = \code{NULL}. If the covariates
#' contain mixed data types, please make sure all categorical variables are properly formatted as \code{factor} type.
#' @param COV.MY a \code{data.frame} or \code{matrix} of covariates dataset for testing the association \code{Y ~ M}.
#' Covariates specified here will not participate penalization. If not specified, the same set of covariates for
#' \code{M ~ X} will be applied (i.e., \code{COV.XM}. Using different sets of covariates is allowed but this needs to be 
#' handled carefully.
#' @param Y.type data type of outcome (\code{Y}). Either \code{'continuous'} (default) or \code{'binary'}.
#' @param M.type data type of mediator (\code{M}). Either \code{'gaussian'} (default) or \code{'negbin'} (i.e., negative binomial).
#' @param penalty the penalty to be applied to the model. Either \code{'MCP'} (the default), \code{'SCAD'}, or
#' \code{'lasso'}.
#' @param topN an integer specifying the number of top markers from sure independent screening.
#' Default = \code{NULL}. If \code{NULL}, \code{topN} will be either \code{ceiling(n/log(n))} for continuous outcome,
#' or \code{ceiling(n/(2*log(n)))} for binary outcome, where \code{n} is the sample size. If the sample size is greater
#' than topN (pre-specified or calculated), all mediators will be included in the test (i.e. low-dimensional scenario).
#' @param scale logical. Should the function scale the data? Default = \code{TRUE}.
#' @param Bonfcut Bonferroni-corrected p value cutoff applied to select significant mediators. Default = \code{0.05}.
#' @param verbose logical. Should the function be verbose? Default = \code{FALSE}.
#' @param parallel logical. Enable parallel computing feature? Default = \code{FALSE}.
#' @param ncore number of cores to run parallel computing Valid when \code{parallel = TRUE}.
#' @param ... other arguments passed to \code{ncvreg}.
#'
#' @return A data.frame containing mediation testing results of selected mediators.
#' \describe{
#'     \item{Index: }{mediation name of selected significant mediator.}
#'     \item{alpha_hat: }{coefficient estimates of exposure (X) --> mediators (M) (adjusted for covariates).}
#'     \item{beta_hat: }{coefficient estimates of mediators (M) --> outcome (Y) (adjusted for covariates and exposure).}
#'     \item{IDE: }{mediation (indirect) effect, i.e., alpha*beta.}
#'     \item{rimp: }{relative importance of the mediator.}
#'     \item{pmax: }{joint raw p-value of selected significant mediator (based on Bonferroni method).}
#' }
#'
#' @references Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J, Just A, Colicino E, Vokonas P, Zhao L,
#' Lv J, Baccarelli A, Hou L, Liu L. Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies.
#' Bioinformatics. 2016. DOI: 10.1093/bioinformatics/btw351. PMID: 27357171; PMCID: PMC5048064
#'
#' @examples
#' \dontrun{
#' # Note: In the following examples, M1, M2, and M3 are true mediators.

#'
#' # When Y is continuous and normally distributed
#' # Example 1 (continuous outcome):
#' data(ContinuousOutcome)
#' pheno_data <- ContinuousOutcome$PhenoData
#' mediator_data <- ContinuousOutcome$Mediator
#'
#' hima.fit <- hima_classic(
#'   X = pheno_data$Treatment,
#'   Y = pheno_data$Outcome,
#'   M = mediator_data,
#'   COV.XM = pheno_data[, c("Sex", "Age")],
#'   Y.type = "continuous",
#'   scale = FALSE, # Disabled only for simulation data
#'   verbose = TRUE
#' )
#' hima.fit
#'
#' # When Y is binary
#' # Example 2 (binary outcome):
#' data(BinaryOutcome)
#' pheno_data <- BinaryOutcome$PhenoData
#' mediator_data <- BinaryOutcome$Mediator
#'
#' hima.logistic.fit <- hima_classic(
#'   X = pheno_data$Treatment,
#'   Y = pheno_data$Disease,
#'   M = mediator_data,
#'   COV.XM = pheno_data[, c("Sex", "Age")],
#'   Y.type = "binary",
#'   scale = FALSE, # Disabled only for simulation data
#'   verbose = TRUE
#' )
#' hima.logistic.fit
#' }
#'
#' @export
#'
hima_classic <- function(X, M, Y, COV.XM = NULL, COV.MY = COV.XM,
                        Y.type = c("continuous", "binary"),
                        M.type = c("gaussian", "negbin"),
                        penalty = c("MCP", "SCAD", "lasso"),
                        topN = NULL,
                        scale = TRUE,
                        Bonfcut = 0.05,
                        verbose = FALSE,
                        parallel = FALSE,
                        ncore = 1,
                        ...) {
  
  if (!is.null(COV.XM) && !is.data.frame(COV.XM) && !is.matrix(COV.XM)) {
    stop("COV.XM must be either a data.frame or a matrix when provided.")
  }
  if (!is.null(COV.MY) && !is.data.frame(COV.MY) && !is.matrix(COV.MY)) {
    stop("COV.MY must be either a data.frame or a matrix when provided.")
  }
  
  Y.type <- match.arg(Y.type)
  Y.family <- switch(Y.type,
    continuous = "gaussian",
    binary = "binomial",
    stop("Invalid Y.type. Expected 'continuous' or 'binary'.")
  )

  M.type <- match.arg(M.type)
  penalty <- match.arg(penalty)
  
  if (Y.family == "gaussian" && verbose) message("Running linear HIMA with ", penalty, " penalty...")
  if (Y.family == "binomial" && verbose) message("Running logistic HIMA with ", penalty, " penalty...")
  
  if (parallel && (ncore == 1)) ncore <- parallel::detectCores()
  if (!parallel && (ncore > 1)) parallel <- TRUE

  n <- nrow(M)
  p <- ncol(M)

  # Process required variables
  X <- process_var(X, scale)
  M <- process_var(M, scale)

  # Process optional covariates
  COV.XM <- process_var(COV.XM, scale)
  COV.MY <- process_var(COV.MY, scale)

  if (scale && verbose) message("Data scaling is completed.")

  if (is.null(topN)) {
    if (Y.type == "binary") d <- ceiling(n / (2 * log(n))) else d <- ceiling(2 * n / log(n))
  } else {
    d <- topN # the number of top mediators that associated with exposure (X)
  }

  d <- min(p, d) # if d > p select all mediators

  #########################################################################
  ################################ STEP 1 #################################
  #########################################################################
  if (verbose) message("Step 1: Sure Independent Screening ...", "     (", format(Sys.time(), "%X"), ")")

  if (Y.type == "binary") {
    # Screen M using X given the limited information provided by Y (binary)
    if (verbose) message("    Screening M using the association between X (independent variable) and \nM (dependent variable): ", appendLF = FALSE)
    alpha <- SIS_Results <- himasis(NA, M, X, COV.XM,
      glm.family = M.type, modelstatement = "Mone ~ X",
      parallel = parallel, ncore = ncore, verbose, tag = paste0("Sure Independent Screening (M ~ X + COV.XM, family: ", M.type, ")")
    )
    SIS_Pvalue <- SIS_Results[2, ]
  } else if (Y.type == "continuous") {
    # Screen M using Y (continuous)
    if (verbose) message("    Screening M using the association between M (independent variable) and \nY (dependent variable): ", appendLF = FALSE)
    SIS_Results <- himasis(Y, M, X, COV.MY,
      glm.family = "gaussian", modelstatement = "Y ~ Mone + X",
      parallel = parallel, ncore = ncore, verbose, tag = paste0("Sure Independent Screening (Y ~ M + X + COV.MY, family: ", Y.type, ")")
    )
    SIS_Pvalue <- SIS_Results[2, ]
  } else {
    stop(paste0("A ", Y.type, " data type for outcome is not supported."))
  }

  # Note: ranking using p on un-standardized data is equivalent to ranking using beta on standardized data
  SIS_Pvalue_sort <- sort(SIS_Pvalue)
  ID <- which(SIS_Pvalue <= SIS_Pvalue_sort[d]) # the index of top mediators
  if (verbose) message("    Top ", length(ID), " mediators are selected: ", paste0(names(SIS_Pvalue_sort[seq_len(d)]), collapse = ", "))

  M_SIS <- M[, ID]
  M_ID_name <- colnames(M)
  XM <- cbind(M_SIS, X)

  #########################################################################
  ################################ STEP 2 #################################
  #########################################################################
  if (verbose) message("Step 2: Penalized estimate (", penalty, ") ...", "     (", format(Sys.time(), "%X"), ")")

  ## Based on the screening results in step 1. We will find the most influential M on Y.
  if (is.null(COV.MY)) {
    fit <- ncvreg(XM, Y,
      family = Y.family,
      penalty = penalty,
      penalty.factor = c(rep(1, ncol(M_SIS)), 0), ...
    )
  } else {
    COV.MY <- data.frame(COV.MY)
    COV.MY <- data.frame(model.matrix(~., COV.MY))[, -1]
    conf.names <- colnames(COV.MY)
    if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    XM_COV <- cbind(XM, COV.MY)
    fit <- ncvreg(XM_COV, Y,
      family = Y.family,
      penalty = penalty,
      penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(as.matrix(COV.MY)))), ...
    )
  }
  lam <- fit$lambda[which.min(BIC(fit))]
  if (verbose) message("    Tuning parameter lambda selected: ", lam)
  Coefficients <- coef(fit, lambda = lam)
  est <- Coefficients[2:(d + 1)]
  ID_1_non <- which(est != 0)
  if (length(ID_1_non) == 0) {
    if (verbose) message("    All ", penalty, " beta estimates of the ", length(ID), " mediators are zero.")
    results <- NULL
    return(results)
  } else {
    if (verbose) message("    Non-zero ", penalty, " beta estimate(s) of mediator(s) found: ", paste0(names(ID_1_non), collapse = ","))
    beta_est <- est[ID_1_non] # The non-zero MCP estimators of beta
    ID_test <- ID[ID_1_non] # The index of the ID of non-zero beta in Y ~ M
    ##

    if (Y.type == "binary") {
      ## This has been done in step 1 (when Y is binary, alpha is estimated in M ~ X)
      alpha <- alpha[, ID_test, drop = FALSE]
      if (verbose) message("    Using alpha estimated in Step 1 ...   (", format(Sys.time(), "%X"), ")")
    } else if (Y.type == "continuous") {
      if (verbose) message("    Estimating alpha (effect of X on M): ", appendLF = FALSE)
      alpha <- himasis(NA, M[, ID_test, drop = FALSE], X, COV.XM,
        glm.family = M.type,
        modelstatement = "Mone ~ X", parallel = FALSE, ncore = ncore,
        verbose, tag = paste0("site-by-site alpha estimation (M ~ X + COV.XM, family: ", M.type, ")")
      )
    } else {
      stop(paste0("A ", Y.type, " data type for outcome is not supported."))
    }

    #########################################################################
    ################################ STEP 3 #################################
    #########################################################################
    if (verbose) message("Step 3: Joint significance test ...", "     (", format(Sys.time(), "%X"), ")")

    alpha_est_ID_test <- as.numeric(alpha[1, ]) #  the estimator for alpha
    P_alpha <- alpha[2, ] # the raw p-value for alpha
    alpha_est <- alpha_est_ID_test

    ## Post-test based on the oracle property of the MCP penalty
    if (is.null(COV.MY)) {
      YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X)
    } else {
      YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X, COV.MY)
    }

    if (Y.type == "continuous") res <- summary(glm(Y ~ ., family = "gaussian", data = YMX))$coefficients
    if (Y.type == "binary") res <- summary(glm(Y ~ ., family = "binomial", data = YMX))$coefficients

    est <- res[2:(length(ID_test) + 1), 1] # the estimator for beta
    P_beta <- res[2:(length(ID_test) + 1), 4] # the raw p-value for beta

    IDE <- alpha_est * beta_est # mediation(indirect) effect

    ## Use the maximum value as p value
    Pmax <- apply(cbind(P_alpha, P_beta), 1, max)

    ## Bonferroni
    Pmax_Bonf <- Pmax * length(ID_test)
    sig_ind <- which(Pmax_Bonf < Bonfcut)

    results <- data.frame(
      Index = M_ID_name[ID_test][sig_ind],
      alpha_hat = alpha_est[sig_ind],
      beta_hat = beta_est[sig_ind],
      IDE = IDE[sig_ind],
      rimp = (abs(IDE) / sum(abs(IDE)))[sig_ind] * 100,
      pmax = Pmax[sig_ind], row.names = NULL
    )

    if (verbose) message("Done!", "     (", format(Sys.time(), "%X"), ")")

    doParallel::stopImplicitCluster()

    return(results)
  }
}
