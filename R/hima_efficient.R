# This is the function for efficient high-dimensional mediation analysis
#' Efficient high-dimensional mediation analysis
#'
#' \code{hima_efficient} is used to estimate and test high-dimensional mediation effects using an efficient algorithm. It provides
#' higher statistical power than the standard \code{hima}. Note: efficient HIMA is only applicable to mediators and outcomes that
#' are both continuous and normally distributed.
#'
#' @param X a vector of exposure. Do not use \code{data.frame} or \code{matrix}.
#' @param M a \code{data.frame} or \code{matrix} of high-dimensional mediators. Rows represent samples, columns
#' represent mediator variables. \code{M} has to be continuous and normally distributed.
#' @param Y a vector of continuous outcome. Do not use \code{data.frame} or \code{matrix}.
#' @param COV a matrix of adjusting covariates. Rows represent samples, columns represent variables. Can be \code{NULL}.
#' @param topN an integer specifying the number of top markers from sure independent screening.
#' Default = \code{NULL}. If \code{NULL}, \code{topN} will be \code{2*ceiling(n/log(n))}, where \code{n} is the sample size.
#' If the sample size is greater than topN (pre-specified or calculated), all mediators will be included in the test (i.e. low-dimensional scenario).
#' @param scale logical. Should the function scale the data? Default = \code{TRUE}.
#' @param FDRcut Benjamini-Hochberg FDR cutoff applied to select significant mediators. Default = \code{0.05}.
#' @param verbose logical. Should the function be verbose? Default = \code{FALSE}.
#' @param parallel logical. Enable parallel computing feature? Default = \code{FALSE}.
#' @param ncore number of cores to run parallel computing Valid when \code{parallel = TRUE}.
#'
#' @return A data.frame containing mediation testing results of significant mediators (FDR <\code{FDRcut}).
#' \describe{
#'     \item{Index: }{mediation name of selected significant mediator.}
#'     \item{alpha_hat: }{coefficient estimates of exposure (X) --> mediators (M) (adjusted for covariates).}
#'     \item{alpha_se: }{standard error for alpha.}
#'     \item{beta_hat: }{coefficient estimates of mediators (M) --> outcome (Y) (adjusted for covariates and exposure).}
#'     \item{beta_se: }{standard error for beta.}
#'     \item{IDE: }{mediation (indirect) effect, i.e., alpha*beta.}
#'     \item{rimp: }{relative importance of the mediator.}
#'     \item{pmax: }{joint raw p-value of selected significant mediator (based on divide-aggregate composite-null test [DACT] method).}
#' }
#'
#' @references Bai X, Zheng Y, Hou L, Zheng C, Liu L, Zhang H. An Efficient Testing Procedure for High-dimensional Mediators with FDR Control.
#' Statistics in Biosciences. 2024. DOI: 10.1007/s12561-024-09447-4.
#'
#' @examples
#' \dontrun{
#' # Note: In the following example, M1, M2, and M3 are true mediators.
#'
#' # Y is continuous and normally distributed
#' # Example (continuous outcome):
#' data(ContinuousOutcome)
#' pheno_data <- ContinuousOutcome$PhenoData
#' mediator_data <- ContinuousOutcome$Mediator
#'
#' hima_efficient.fit <- hima_efficient(
#'   X = pheno_data$Treatment,
#'   Y = pheno_data$Outcome,
#'   M = mediator_data,
#'   COV = pheno_data[, c("Sex", "Age")],
#'   scale = FALSE, # Disabled only for simulation data
#'   FDRcut = 0.05,
#'   verbose = TRUE
#' )
#' hima_efficient.fit
#' }
#'
#' @export
hima_efficient <- function(X, M, Y, COV = NULL,
                  topN = NULL,
                  scale = TRUE,
                  FDRcut = 0.05,
                  verbose = FALSE,
                  parallel = FALSE,
                  ncore = 1) {
  n <- nrow(M)
  p <- ncol(M)

  # Process required variables
  X <- process_var(X, scale)
  M <- process_var(M, scale)

  # Process optional covariates
  COV <- process_var(COV, scale)

  if (scale && verbose) message("Data scaling is completed.")

  checkParallel("hima_efficient", parallel, ncore, verbose)

  if (is.null(COV)) {
    MZX <- cbind(M, X)
    XZ <- X
    q <- 0
  } else {
    MZX <- cbind(M, COV, X)
    XZ <- cbind(X, COV)
    q <- ncol(COV)
  }

  if (is.null(topN)) d <- ceiling(2 * n / log(n)) else d <- topN # the number of top mediators that associated with exposure (X)
  d <- min(p, d) # if d > p select all mediators

  M_ID_name <- colnames(M)
  if (is.null(M_ID_name)) M_ID_name <- seq_len(p)

  #------------- Step 1:  mediator screening ---------------------------
  if (verbose) message("Step 1: Sure Independent Screening + minimax concave penalty (MCP) ...", "     (", format(Sys.time(), "%X"), ")")

  beta_SIS <- foreach(i = seq_len(p), .combine = "c") %dopar% {
    ID_S <- c(i, (p + 1):(p + q + 1))
    MZX_SIS <- MZX[, ID_S]
    fit <- lsfit(MZX_SIS, Y, intercept = TRUE)
    fit$coefficients[2]
  }

  ## est_a for SIS #########
  alpha_SIS <- foreach(i = seq_len(p), .combine = "c") %dopar% {
    fit_a <- lsfit(XZ, M[, i], intercept = TRUE)
    as.numeric(coef(fit_a))[2]
  }
  ab_SIS <- alpha_SIS * beta_SIS
  ID_SIS <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[d]) # \Omega_1
  d <- length(ID_SIS)

  if (verbose) message("        Top ", d, " mediators are selected: ", paste0(M_ID_name[ID_SIS], collapse = ", "))

  if (verbose) {
    if (is.null(COV)) {
      message("        No covariate was adjusted.")
    } else {
      message("        Adjusting for covariate(s): ", paste0(colnames(COV), collapse = ", "))
    }
  }

  MZX_SIS <- MZX[, c(ID_SIS, (p + 1):(p + q + 1))] # select m_i in \Omega_1 from M
  fit <- ncvreg(MZX_SIS, Y, family = "gaussian", penalty = "MCP")
  lam <- fit$lambda[which.min(BIC(fit))]
  beta_penalty <- coef(fit, lambda = lam)[2:(d + 1)]
  id_non <- ID_SIS[which(beta_penalty != 0)] # the ID of non-zero

  #----------- Step 2: Refitted partial regression ----------------------
  if (verbose) message("Step 2: Refitted partial regression ...", "     (", format(Sys.time(), "%X"), ")")
  ## beta_est ########
  MZX_penalty <- MZX[, c(id_non, (p + 1):(p + q + 1))]
  fit <- lsfit(MZX_penalty, Y, intercept = TRUE)
  beta_est_cox <- fit$coefficients[2:(dim(MZX_penalty)[2] + 1)]
  beta_SE_cox <- ls.diag(fit)$std.err[2:(dim(MZX_penalty)[2] + 1)]

  # Computes basic statistics, including standard errors, t- and p-values for the regression coefficients.
  beta_est <- fit$coefficients[2:(length(id_non) + 1)] # estimated beta != 0
  beta_SE <- ls.diag(fit)$std.err[2:(length(id_non) + 1)]
  P_beta_penalty <- 2 * (1 - pnorm(abs(beta_est_cox[seq_along(id_non)]) / beta_SE_cox[seq_along(id_non)], 0, 1))

  beta_results <- foreach(i = seq_len(p), .combine = rbind) %dopar% {
    if (i %in% id_non) {
      j <- which(id_non == i)
      c(beta_est[j], beta_SE[j], P_beta_penalty[j])
    } else {
      MZX_ora <- MZX[, c(id_non, i, (p + 1):(p + q + 1))]
      fit_ora <- lsfit(MZX_ora, Y, intercept = TRUE)
      beta_est_cox <- fit_ora$coefficients[2:(ncol(MZX_ora) + 1)]
      est <- beta_est_cox[length(id_non) + 1]
      beta_se_cox <- ls.diag(fit_ora)$std.err[2:(ncol(MZX_ora) + 1)]
      se <- beta_se_cox[length(id_non) + 1]
      p_val <- 2 * (1 - pnorm(abs(beta_est_cox) / beta_se_cox, 0, 1))[length(id_non) + 1]
      c(est, se, p_val)
    }
  }
  if (is.null(dim(beta_results))) beta_results <- matrix(beta_results, nrow = 1)
  beta_est_orc <- beta_results[, 1]
  beta_SE_orc <- beta_results[, 2]
  P_oracle_beta <- beta_results[, 3]

  ## ----- P_oracle_alpha ----------------- ##
  alpha_pen_results <- foreach(i = seq_along(id_non), .combine = rbind) %dopar% {
    fit_a <- lsfit(XZ, M[, id_non[i]], intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    se_a <- ls.diag(fit_a)$std.err[2]
    sd_1 <- abs(est_a) / se_a
    c(est_a, se_a, 2 * (1 - pnorm(sd_1, 0, 1)))
  }
  if (is.null(dim(alpha_pen_results))) alpha_pen_results <- matrix(alpha_pen_results, nrow = 1)
  alpha_est_penalty <- alpha_pen_results[, 1]
  alpha_SE_penalty <- alpha_pen_results[, 2]
  P_alpha_penalty <- alpha_pen_results[, 3]

  ### P_oracle_alpha #########
  P_oracle_alpha <- matrix(0, 1, p) # an empty vector
  alpha_est_orc <- matrix(0, 1, p)
  alpha_SE_orc <- matrix(0, 1, p)

  alpha_results <- foreach(i = seq_len(p), .combine = rbind) %dopar% {
    if (i %in% id_non) {
      j <- which(id_non == i)
      c(alpha_est_penalty[j], alpha_SE_penalty[j], P_alpha_penalty[j])
    } else {
      fit_a_ora <- lsfit(XZ, M[, i], intercept = TRUE)
      est_a_ora <- matrix(coef(fit_a_ora))[2]
      se_a_ora <- ls.diag(fit_a_ora)$std.err[2]
      sd_1_ora <- abs(est_a_ora) / se_a_ora
      p_val <- 2 * (1 - pnorm(sd_1_ora, 0, 1))
      c(est_a_ora, se_a_ora, p_val)
    }
  }
  if (is.null(dim(alpha_results))) alpha_results <- matrix(alpha_results, nrow = 1)
  alpha_est_orc <- alpha_results[, 1]
  alpha_SE_orc <- alpha_results[, 2]
  P_oracle_alpha <- alpha_results[, 3]

  #---------- Step 3: DACT  -------------------------
  if (verbose) message("Step 3: Divide-aggregate composite-null test (DACT) ...", "     (", format(Sys.time(), "%X"), ")")

  # Mediator selection
  P_oracle_alpha[P_oracle_alpha == 0] <- 10^(-17)
  P_oracle_beta[P_oracle_beta == 0] <- 10^(-17)
  P_BH <- (1:p) * (FDRcut / p)

  ## DACT
  DACT_ora <- DACT(p_a = t(P_oracle_alpha), p_b = t(P_oracle_beta))
  P_sort_DACT <- sort(DACT_ora)
  SN <- sum(as.numeric(P_sort_DACT <= P_BH))
  ID_BH_DACT <- order(DACT_ora)[1:SN]

  # # Total effect
  # if(is.null(COV)) {
  #   YX <- data.frame(Y = Y, X = X)
  # } else {
  #   YX <- data.frame(Y = Y, X = X, COV)
  # }
  #
  # gamma_est <- coef(glm(Y ~ ., family = Y.family, data = YX))[2]

  IDE <- alpha_est_orc[ID_BH_DACT] * beta_est_orc[ID_BH_DACT]

  if (length(ID_BH_DACT) > 0) {
    out_result <- data.frame(
      Index = M_ID_name[ID_BH_DACT],
      alpha_hat = alpha_est_orc[ID_BH_DACT],
      alpha_se = alpha_SE_orc[ID_BH_DACT],
      beta_hat = beta_est_orc[ID_BH_DACT],
      beta_se = beta_SE_orc[ID_BH_DACT],
      IDE = IDE,
      rimp = abs(IDE) / sum(abs(IDE)) * 100,
      pmax = DACT_ora[ID_BH_DACT], row.names = NULL
    )
    if (verbose) message(paste0("        ", length(ID_BH_DACT), " significant mediator(s) identified."))
  } else {
    if (verbose) message("        No significant mediator identified.")
    out_result <- NULL
  }

  if (verbose) message("Done!", "     (", format(Sys.time(), "%X"), ")")

  doParallel::stopImplicitCluster()
  return(out_result)
}
