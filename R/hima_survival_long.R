# This is the function for our proposed method for high-dimensional and longitudinal mediation analysis for survival data.
#' High-dimensional mediation analysis for longitudinal mediator and survival outcome data
#'
#' \code{hima_survival_long} estimates and tests high-dimensional longitudinal mediation effects for survival data in a counting
#' process framework.
#'
#' @param X A numeric vector of exposure values (do not use \code{data.frame} or \code{matrix}).
#' @param M A \code{data.frame} or \code{matrix} of high-dimensional mediators (rows = observations/intervals, columns = mediators).
#' @param tstart A numeric vector of starting times for each observation/interval (e.g., entry time in a counting-process setup).
#' @param tstop A numeric vector of stopping times for each observation/interval (e.g., event/censoring time in a counting-process setup).
#' @param status A numeric vector of censoring indicators (\code{1 = event}, \code{0 = censored}).
#' @param id A vector of subject identifiers (used for clustering/random effects).
#' @param COV A \code{matrix} or \code{data.frame} of adjusting covariates. Rows represent samples, columns represent variables. Can be \code{NULL}.
#' @param topN Integer specifying the number of top mediators retained after sure independent screening (SIS). If \code{NULL}
#' (default), \code{topN = ceiling(n/log(n))}, where \code{n} is the number of unique subjects. When \code{topN} exceeds the
#' total number of mediators, all mediators are kept (i.e., the low-dimensional scenario).
#' @param scale Logical. Should the function scale the exposure, mediators, and covariates? Default = \code{TRUE}.
#' @param Bonfcut Bonferroni-corrected p value cutoff applied to select significant mediators. Default = \code{0.05}.
#' @param verbose Logical. Should progress messages be printed? Default = \code{FALSE}.
#' @param parallel Logical. Enable parallel computing for SIS? Default = \code{FALSE}.
#' @param ncore Integer specifying the number of cores to use when \code{parallel = TRUE}.
#'
#' @return A \code{data.frame} containing mediation testing results of significant mediators (joint p-value < \code{Bonfcut}).
#' \describe{
#'   \item{Index}{Mediator name of the selected significant mediator.}
#'   \item{alpha_hat}{Coefficient estimates for the exposure (X) --> mediator (M) model (adjusted for covariates).}
#'   \item{alpha_se}{Standard error for \code{alpha_hat}.}
#'   \item{beta_hat}{Coefficient estimates for the mediator (M) --> outcome (Y) model (adjusted for covariates and exposure).}
#'   \item{beta_se}{Standard error for \code{beta_hat}.}
#'   \item{IDE}{Indirect (mediation) effect estimate, i.e., \code{alpha_hat * beta_hat}.}
#'   \item{rimp}{Relative importance of the mediator.}
#'   \item{pmax}{joint raw p-value of selected significant mediator (based on Bonferroni method).}
#' }
#'
#' @references
#' Liu L, Zhang H, Zheng Y, Gao T, Zheng C, Zhang K, Hou L, Liu L. High-dimensional mediation analysis for longitudinal mediators and 
#' survival outcomes. Briefings in Bioinformatics. 2025. DOI: 10.1093/bib/bbaf206. PMID: 40350699 PMCID: PMC12066418
#'
#' @examples
#' \dontrun{
#' data(SurvivalLongData)
#' pheno_data <- SurvivalLongData$PhenoData
#' mediator_data <- SurvivalLongData$Mediator
#'
#' hima_survival_long.fit <- hima_survival_long(
#'   X = pheno_data$Treatment,
#'   M = mediator_data,
#'   tstart = pheno_data$Tstart,
#'   tstop = pheno_data$Tstop,
#'   status = pheno_data$Status,
#'   id = pheno_data$ID,
#'   COV = pheno_data[, c("Sex", "Age")],
#'   verbose = TRUE
#' )
#' hima_survival_long.fit
#' }
#'
#' @export
hima_survival_long <- function(X,
                               M,
                               tstart,
                               tstop,
                               status,
                               id,
                               COV = NULL,
                               topN = NULL,
                               scale = TRUE,
                               Bonfcut = 0.05,
                               verbose = FALSE,
                               parallel = FALSE,
                               ncore = 1) {
  X <- as.numeric(X)
  M <- as.matrix(M)
  tstart <- as.numeric(tstart)
  tstop <- as.numeric(tstop)
  status <- as.numeric(status)

  if (any(lengths(list(X, tstart, tstop, status, id)) != nrow(M))) {
    stop("X, tstart, tstop, status, id, and mediators must have the same number of rows/observations.")
  }
  
  if (!is.numeric(Bonfcut) || Bonfcut <= 0) {
    stop("'Bonfcut' must be a positive numeric value.")
  }
  Bonfcut <- min(Bonfcut, 1)
  
  mediator_labels <- colnames(M)
  if (is.null(mediator_labels)) mediator_labels <- as.character(seq_len(ncol(M)))
  mediator_model_names <- make.names(mediator_labels, unique = TRUE)
  
  if (is.null(COV)) {
    q <- 0
    covariate_names <- character(0)
  } else {
    COV <- as.matrix(COV)
    if (nrow(COV) != nrow(M)) stop("COV must have the same number of rows as M.")
    if (is.null(colnames(COV))) {
      colnames(COV) <- paste0("COV", seq_len(ncol(COV)))
    }
    covariate_names <- colnames(COV)
    q <- ncol(COV)
  }
  
  M <- process_var(M, scale)
  if (!is.null(COV)) {
    COV <- process_var(COV, scale)
  }
  X <- as.numeric(process_var(X, scale))
  
  if (scale && verbose) message("Data scaling is completed.")
  checkParallel("hima_survival_long", parallel, ncore, verbose)
  
  colnames(M) <- mediator_model_names
  if (!is.null(COV)) {
    colnames(COV) <- covariate_names
  }
  
  n_subject <- length(unique(id))
  p <- ncol(M)
  
  if (is.null(topN)) {
    divisor <- max(log(max(n_subject, 2)), 1)
    d_0 <- ceiling(n_subject / (2 * divisor))
  } else {
    d_0 <- as.integer(topN)
  }
  d_0 <- max(1, min(p, d_0))
  
  if (verbose) {
    message("Step 1: Sure Independent Screening ...     (", format(Sys.time(), "%X"), ")")
  }
  
  base_cox <- data.frame(
    tstart = tstart,
    tstop = tstop,
    status = status,
    X = X,
    id = id,
    stringsAsFactors = FALSE
  )
  base_lme <- data.frame(
    id = id,
    X = X,
    stringsAsFactors = FALSE
  )
  if (!is.null(COV)) {
    cov_df <- as.data.frame(COV)
    base_cox <- cbind(base_cox, cov_df)
    base_lme <- cbind(base_lme, cov_df)
  }
  
  quote_terms <- function(x) {
    if (length(x) == 0) {
      return(character(0))
    }
    needs_quote <- make.names(x) != x
    x[needs_quote] <- paste0("`", x[needs_quote], "`")
    x
  }
  
  cox_terms <- c("X", covariate_names, "mediator")
  cox_formula <- stats::as.formula(
    paste("survival::Surv(tstart, tstop, status) ~", paste(quote_terms(cox_terms), collapse = " + "))
  )
  
  lme_terms <- c("X", covariate_names)
  lme_formula <- stats::as.formula(
    paste("mediator ~", paste(quote_terms(lme_terms), collapse = " + "))
  )
  
  sis_results <- foreach::foreach(
    k = seq_len(p),
    .combine = rbind,
    .packages = c("survival", "nlme")
  ) %dopar% {
    Mk <- M[, k]
    cox_data <- cbind(base_cox, mediator = Mk)
    beta_est <- tryCatch({
      fit1 <- survival::coxph(cox_formula, data = cox_data, cluster = cox_data$id)
      stats::coef(fit1)["mediator"]
    }, error = function(e) NA_real_)
    if (is.na(beta_est)) beta_est <- 0
    
    lme_data <- cbind(base_lme, mediator = Mk)
    lme_fit <- tryCatch({
      nlme::lme(
        fixed = lme_formula,
        random = ~1 | id,
        data = lme_data,
        control = nlme::lmeControl(opt = "optim")
      )
    }, error = function(e) NULL)
    
    if (is.null(lme_fit)) {
      alpha_est <- 0
      alpha_se <- 0
      alpha_p <- 1
    } else {
      t_table <- summary(lme_fit)$tTable
      alpha_est <- t_table["X", "Value"]
      alpha_se <- t_table["X", "Std.Error"]
      alpha_p <- t_table["X", "p-value"]
    }
    
    c(beta_est, alpha_est, alpha_se, alpha_p)
  }
  colnames(sis_results) <- c("beta", "alpha", "alpha_se", "alpha_p")
  
  ab_SIS <- sis_results[, "beta"] * sis_results[, "alpha"]
  ab_SIS[is.na(ab_SIS)] <- 0
  ID_SIS <- order(-abs(ab_SIS))[seq_len(d_0)]
  
  if (verbose) {
    message("        Top ", length(ID_SIS), " mediators are selected: ", paste0(mediator_labels[ID_SIS], collapse = ", "))
    message("Step 2: L1-penalized Cox lasso estimates ...     (", format(Sys.time(), "%X"), ")")
  }
  
  design_mat <- as.matrix(cbind(X, COV, M[, ID_SIS, drop = FALSE]))
  penalty_factors <- c(rep(0, 1 + q), rep(1, length(ID_SIS)))
  
  cv_fit <- glmnet::cv.glmnet(
    x = design_mat,
    y = survival::Surv(tstart, tstop, status),
    family = "cox",
    penalty.factor = penalty_factors
  )
  fit_lasso <- glmnet::glmnet(
    x = design_mat,
    y = survival::Surv(tstart, tstop, status),
    family = "cox",
    lambda = cv_fit$lambda.min,
    penalty.factor = penalty_factors
  )
  lasso_coef <- as.vector(stats::coef(fit_lasso))[-seq_len(q + 1)]
  lasso_coef[is.na(lasso_coef)] <- 0
  ID_lasso <- ID_SIS[which(lasso_coef != 0)]
  
  if (length(ID_lasso) == 0) {
    if (verbose) message("        No mediator passed the Lasso selection step.")
    return(NULL)
  }
  
  if (verbose) {
    message("        Top ", length(ID_lasso), " mediators are selected: ", paste0(mediator_labels[ID_lasso], collapse = ", "))
    message("Step 3: Multiple-testing procedure ...     (", format(Sys.time(), "%X"), ")")
  }
  
  cox_lasso_data <- cbind(base_cox, M[, ID_lasso, drop = FALSE])
  cox_terms_full <- c("X", covariate_names, mediator_model_names[ID_lasso])
  cox_formula_full <- stats::as.formula(
    paste("survival::Surv(tstart, tstop, status) ~", paste(quote_terms(cox_terms_full), collapse = " + "))
  )
  fit_full <- survival::coxph(cox_formula_full, data = cox_lasso_data, cluster = cox_lasso_data$id)
  coef_summary <- summary(fit_full)$coef
  beta_block <- coef_summary[mediator_model_names[ID_lasso], , drop = FALSE]
  beta_hat <- beta_block[, "coef"]
  beta_se <- beta_block[, "se(coef)"]
  beta_pvalue <- beta_block[, "Pr(>|z|)"]
  
  bind_p <- rbind(
    pmin(beta_pvalue * length(ID_lasso), 1),
    pmin(sis_results[ID_lasso, "alpha_p"] * length(ID_lasso), 1)
  )
  final_p <- apply(bind_p, 2, max)
  selected <- which(final_p <= Bonfcut)
  
  if (length(selected) == 0) {
    if (verbose) message("        No significant mediator identified.")
    return(NULL)
  }
  
  ID_fdr <- ID_lasso[selected]
  
  beta_hat <- beta_hat[selected]
  beta_se <- beta_se[selected]
  alpha_hat <- sis_results[ID_fdr, "alpha"]
  alpha_se <- sis_results[ID_fdr, "alpha_se"]
  IDE <- alpha_hat * beta_hat
  
  if (sum(abs(IDE)) == 0) {
    rimp <- rep(0, length(IDE))
  } else {
    rimp <- abs(IDE) / sum(abs(IDE)) * 100
  }
  
  out_result <- data.frame(
    Index = mediator_labels[ID_fdr],
    alpha_hat = alpha_hat,
    alpha_se = alpha_se,
    beta_hat = beta_hat,
    beta_se = beta_se,
    IDE = IDE,
    rimp = rimp,
    pmax = final_p[selected],
    row.names = NULL
  )
  
  if (verbose) {
    message(
      paste0(
        "        ", length(ID_fdr), " significant mediator(s) identified: ",
        paste0(mediator_labels[ID_fdr], collapse = ", "
        )
      )
    )
    message("Done!     (", format(Sys.time(), "%X"), ")")
  }
  
  return(out_result)
}
