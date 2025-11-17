# This is the function for high-dimensional quantile mediation analysis
#' High-dimensional quantile mediation analysis
#'
#' \code{hima_quantile} is used to estimate and test high-dimensional quantile mediation effects.
#'
#' @param X a vector of exposure. Do not use \code{data.frame} or \code{matrix}.
#' @param M a \code{data.frame} or \code{matrix} of high-dimensional mediators. Rows represent samples, columns
#' represent mediator variables.
#' @param Y a vector of continuous outcome. Do not use \code{data.frame} or \code{matrix}.
#' @param COV a matrix of adjusting covariates. Rows represent samples, columns represent variables. Can be \code{NULL}.
#' @param penalty the penalty to be applied to the model (a parameter passed to function \code{conquer.cv.reg} in package \code{conquer}.
#' Either \code{'MCP'} (the default), \code{'SCAD'}, or \code{'lasso'}.
#' @param topN an integer specifying the number of top markers from sure independent screening.
#' Default = \code{NULL}. If \code{NULL}, \code{topN} will be \code{2*ceiling(n/log(n))}, where \code{n} is the sample size.
#' If the sample size is greater than topN (pre-specified or calculated), all mediators will be included in the test (i.e. 
#' low-dimensional scenario).
#' @param tau quantile level of outcome. Default = \code{0.5}. A vector of tau is accepted.
#' @param scale logical. Should the function scale the data? Default = \code{TRUE}.
#' @param Bonfcut Bonferroni-corrected p value cutoff applied to select significant mediators. Default = \code{0.05}.
#' @param verbose logical. Should the function be verbose? Default = \code{FALSE}.
#' @param parallel logical. Enable parallel computing feature? Default = \code{FALSE}.
#' @param ncore number of cores to run parallel computing Valid when \code{parallel = TRUE}.
#' @param ... reserved passing parameter.
#'
#' @return A data.frame containing mediation testing results of selected mediators (Bonferroni-adjusted p value <\code{Bonfcut}).
#' \describe{
#'     \item{Index: }{mediation name of selected significant mediator.}
#'     \item{alpha_hat: }{coefficient estimates of exposure (X) --> mediators (M) (adjusted for covariates).}
#'     \item{alpha_se: }{standard error for alpha.}
#'     \item{beta_hat: }{coefficient estimates of mediators (M) --> outcome (Y) (adjusted for covariates and exposure).}
#'     \item{beta_se: }{standard error for beta.}
#'     \item{IDE: }{mediation (indirect) effect, i.e., alpha*beta.}
#'     \item{rimp: }{relative importance of the mediator.}
#'     \item{pmax: }{joint raw p-value of selected significant mediator (based on Bonferroni method).}
#' }
#'
#' @references Zhang H, Hong X, Zheng Y, Hou L, Zheng C, Wang X, Liu L. High-Dimensional Quantile Mediation Analysis with Application to a Birth
#' Cohort Study of Mother–Newborn Pairs. Bioinformatics. 2024. DOI: 10.1093/bioinformatics/btae055. PMID: 38290773; PMCID: PMC10873903
#'
#' @examples
#' \dontrun{
#' # Note: In the following example, M1, M2, and M3 are true mediators.
#'
#' data(QuantileData)
#' pheno_data <- QuantileData$PhenoData
#' mediator_data <- QuantileData$Mediator
#'
#' hima_quantile.fit <- hima_quantile(
#'   X = pheno_data$Treatment,
#'   Y = pheno_data$Outcome,
#'   M = mediator_data,
#'   COV = pheno_data[, c("Sex", "Age")],
#'   tau = c(0.3, 0.5, 0.7),
#'   scale = FALSE, # Disabled only for simulation data
#'   Bonfcut = 0.05,
#'   verbose = TRUE
#' )
#' hima_quantile.fit
#' }
#'
#' @export
hima_quantile <- function(X, M, Y, COV = NULL,
                  penalty = c("MCP", "SCAD", "lasso"),
                  topN = NULL,
                  tau = 0.5,
                  scale = TRUE,
                  Bonfcut = 0.05,
                  verbose = FALSE,
                  parallel = FALSE,
                  ncore = 1,
                  ...) {
  penalty <- match.arg(penalty)

  ## --- Validate tau: must be finite numeric in (0, 1) ---
  if (!is.numeric(tau) || any(!is.finite(tau))) {
    stop("`tau` must be a finite numeric vector.")
  }
  if (any(tau <= 0 | tau >= 1)) {
    stop("All `tau` values must lie strictly between 0 and 1.")
  }
  
  ## Outcome checks: only continuous Y is allowed for now
  if (is.factor(Y) && nlevels(Y) == 2L ||
      is.numeric(Y) && all(stats::na.omit(Y) %in% c(0, 1))) {
    stop("hima_quantile() only supports continuous outcomes. ")
  }
  
  if (!is.numeric(Y)) {
    stop("Y must be numeric for hima_quantile().")
  }
  
  n <- nrow(M)
  p <- ncol(M)

  # Process required variables
  X <- process_var(X, scale)
  M <- process_var(M, scale)

  # Process optional covariates
  COV <- process_var(COV, scale)

  if (nrow(X) != n || length(Y) != n || (!is.null(COV) && nrow(COV) != n)) {
    stop("X, M, Y, and COV (if provided) must have the same number of observations.")
  }
  
  if (scale && verbose) message("Data scaling is completed.")

  checkParallel("hima_quantile", parallel, ncore, verbose)

  if (is.null(COV)) XZ <- X else XZ <- cbind(X, COV)

  if (is.null(topN)) d <- ceiling(2 * n / log(n)) else d <- topN # the number of top mediators that associated with exposure (X)
  d <- min(p, d) # if d > p select all mediators

  M_ID_name <- colnames(M)
  if (is.null(M_ID_name)) M_ID_name <- seq_len(p)

  out_result <- NULL

  for (tau_temp in tau)
  {
    if (verbose) message("Running penalized quantile regression with tau = ", tau_temp, " ...", "     (", format(Sys.time(), "%X"), ")")

    #------------- Step 1: Mediator screening ---------------------------
    if (verbose) message("Step 1: Sure Independent Screening ...", "     (", format(Sys.time(), "%X"), ")")

    alpha_est <- matrix(0, 1, p) # the OLS estimator of alpha
    alpha_SE <- matrix(0, 1, p) # the SE of alpha-OLS
    beta_SIS_est <- matrix(0, 1, p) # the screening based estimator of beta
    # beta_SIS_SE <- matrix(0, 1, p) # the SE of beta_SIS_est

    screening_res <- foreach(k = seq_len(p), .combine = rbind) %dopar% {
      MXZ_k <- cbind(M[, k], XZ)
      fit_rq <- rq(Y ~ MXZ_k, tau = tau_temp, method = "fn", model = TRUE)
      beta_val <- fit_rq$coefficients[2]
      fit_M <- lsfit(XZ, M[, k], intercept = TRUE)
      alpha_val <- matrix(coef(fit_M))[2]
      alpha_se <- ls.diag(fit_M)$std.err[2]
      c(beta_val, alpha_val, alpha_se)
    }
    if (is.null(dim(screening_res))) screening_res <- matrix(screening_res, nrow = 1)
    beta_SIS_est <- screening_res[, 1]
    alpha_est <- screening_res[, 2]
    alpha_SE <- screening_res[, 3]

    T_sobel <- beta_SIS_est
    ID_SIS <- which(-abs(T_sobel) <= sort(-abs(T_sobel))[d]) # the index set in Step 1 after the screening

    if (verbose) message("        Top ", length(ID_SIS), " mediators are selected: ", paste0(M_ID_name[ID_SIS], collapse = ", "))


    #----------- Step 2: Penalized estimate in Quantile Regression model
    if (verbose) message("Step 2: Penalized estimate (", penalty, ") ...", "     (", format(Sys.time(), "%X"), ")")

    if (verbose) {
      if (is.null(COV)) {
        message("        No covariate was adjusted.")
      } else {
        message("        Adjusting for covariate(s): ", paste0(colnames(COV), collapse = ", "))
      }
    }

    # Design matrix including only SIS-selected mediators
    MXZ_ID_SIS <- cbind(M[, ID_SIS], XZ)

    fit.penalty <- conquer::conquer.cv.reg(X = MXZ_ID_SIS, Y = Y, tau = tau_temp, penalty = tolower(penalty))
    beta.penalty <- fit.penalty$coeff.min[seq_len(length(ID_SIS)) + 1L]
    
    # Indices of non-zero mediators after penalization
    ID_penalty <- ID_SIS[which(beta.penalty != 0)]
    
    # if no mediators survive penalization, skip Step 3 for this tau
    if (length(ID_penalty) == 0) {
      if (verbose) {
        message("Step 3: Joint significance test ...", "     (", format(Sys.time(), "%X"), ")")
        message("        No mediators selected by penalized regression; skipping joint significance testing for this tau.")
      }
      next
    }
    
    #---------- Step 3: Mediator significance testing
    if (verbose) message("Step 3: Joint significance test ...", "     (", format(Sys.time(), "%X"), ")")

    beta_fit_penalty <- matrix(0, 1, p)
    MXZ_penalty <- cbind(M[, ID_penalty], XZ) # cbind M[ID_penalty], X, and COV

    fit_rq_penalty <- quantreg::rq(Y ~ MXZ_penalty, tau = tau_temp, method = "fn", model = TRUE)
    rq_est_penalty <- summary(fit_rq_penalty, covariance = TRUE, se = "boot", R = 10000)$coefficients
    beta_hat_penalty <- matrix(rq_est_penalty[2:(length(ID_penalty) + 1), 1])
    beta_SE_penalty <- matrix(rq_est_penalty[2:(length(ID_penalty) + 1), 2])
    beta_fit_penalty[ID_penalty] <- beta_hat_penalty

    #--- the p-values
    P_raw_2k <- 2 * (1 - pnorm(abs(alpha_est) / alpha_SE, 0, 1)) # the p-values of alpha
    P_raw_1k_penalty <- 2 * (1 - pnorm(abs(beta_hat_penalty) / beta_SE_penalty, 0, 1))
    P_max_k_penalty <- apply(cbind(P_raw_1k_penalty, P_raw_2k[ID_penalty]), 1, max)

    #-- the Pmax method
    P_penalty_Pmax <- P_max_k_penalty * length(ID_penalty)
    sig_ind <- which(P_penalty_Pmax < Bonfcut)
    ID_Non_penalty_Pmax <- ID_penalty[sig_ind] # the ID of significant M by JS

    IDE <- alpha_est[ID_Non_penalty_Pmax] * beta_hat_penalty[sig_ind]

    if (length(ID_Non_penalty_Pmax) > 0) {
      out_result <- rbind(
        out_result,
        data.frame(
          Index = M_ID_name[ID_Non_penalty_Pmax],
          alpha_hat = alpha_est[ID_Non_penalty_Pmax],
          alpha_se = alpha_SE[ID_Non_penalty_Pmax],
          beta_hat = beta_hat_penalty[sig_ind],
          beta_se = beta_SE_penalty[sig_ind],
          IDE = IDE,
          rimp = abs(IDE) / sum(abs(IDE)) * 100,
          pmax = P_max_k_penalty[sig_ind],
          tau = tau_temp, row.names = NULL
        )
      )
      if (verbose) message(paste0("        ", length(ID_Non_penalty_Pmax), " significant mediator(s) identified."))
    } else {
      if (verbose) message("        No significant mediator identified.")
    }
    if (verbose) message("\t")
  }

  if (verbose) message("Done!", "     (", format(Sys.time(), "%X"), ")")

  doParallel::stopImplicitCluster()
  
  return(out_result)
}
