# This is the function for high-dimensional compositional microbiome mediation analysis
#' High-dimensional mediation analysis for compositional microbiome data
#'
#' \code{hima_microbiome} is used to estimate and test high-dimensional mediation effects for compositional microbiome data.
#'
#' @param X a vector of exposure. Do not use \code{data.frame} or \code{matrix}.
#' @param OTU a \code{data.frame} or \code{matrix} of high-dimensional Operational Taxonomic Unit (OTU) data (mediators).
#' Rows represent samples, columns represent variables.
#' @param Y a vector of continuous outcome. Binary outcome is not allowed. Do not use \code{data.frame} or \code{matrix}.
#' @param COV a \code{data.frame} or \code{matrix} of adjusting covariates. Rows represent samples, columns represent
#' microbiome variables. Can be \code{NULL}.
#' @param FDRcut Hommel FDR cutoff applied to select significant mediators. Default = \code{0.05}.
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
#'     \item{pmax: }{joint raw p-value of selected significant mediator (based on Hommel FDR method).}
#' }
#'
#' @references
#' 1. Zhang H, Chen J, Feng Y, Wang C, Li H, Liu L. Mediation effect selection in high-dimensional and compositional microbiome data.
#' Stat Med. 2021. DOI: 10.1002/sim.8808. PMID: 33205470; PMCID: PMC7855955
#'
#' 2. Zhang H, Chen J, Li Z, Liu L. Testing for mediation effect with application to human microbiome data.
#' Stat Biosci. 2021. DOI: 10.1007/s12561-019-09253-3. PMID: 34093887; PMCID: PMC8177450
#'
#' @examples
#' \dontrun{
#' # Note: In the following example, M1, M2, and M3 are true mediators.
#'
#' data(MicrobiomeData)
#' pheno_data <- MicrobiomeData$PhenoData
#' mediator_data <- MicrobiomeData$Mediator
#' 
#' hima_microbiome.fit <- hima_microbiome(
#'   X = pheno_data$Treatment,
#'   Y = pheno_data$Outcome,
#'   OTU = mediator_data,
#'   COV = pheno_data[, c("Sex", "Age")],
#'   FDRcut = 0.05,
#'   verbose = TRUE
#' )
#' hima_microbiome.fit
#' }
#'
#' @export
hima_microbiome <- function(X,
                            OTU,
                            Y,
                            COV = NULL,
                            FDRcut = 0.05,
                            verbose = FALSE,
                            parallel = FALSE,
                            ncore = 1) {
  
  ## Outcome checks: only continuous Y is allowed for now
  if (is.factor(Y) && nlevels(Y) == 2L ||
      is.numeric(Y) && all(stats::na.omit(Y) %in% c(0, 1))) {
    stop("hima_microbiome() only supports continuous outcomes. ")
  }
  
  if (!is.numeric(Y)) {
    stop("Y must be numeric for hima_microbiome().")
  }
  
  X <- matrix(X, ncol = 1)
  
  M_raw <- as.matrix(OTU)
  
  M_ID_name <- colnames(M_raw)
  if (is.null(M_ID_name)) M_ID_name <- seq_len(ncol(M_raw))
  
  if (!is.null(COV)) {
    COV <- as.matrix(COV)
    X <- cbind(X, COV)
  }
  
  X <- scale(X)
  
  Y <- Y - mean(Y)
  
  M <- M_raw
  n <- dim(M)[1]
  d <- dim(M)[2]
  
  if (nrow(X) != n || length(Y) != n || (!is.null(COV) && nrow(COV) != n)) {
    stop("X, OTU, Y, and COV (if provided) must have the same number of observations.")
  }
  
  alpha_EST <- matrix(0, 1, d)
  alpha_SE <- matrix(0, 1, d)
  beta_EST <- matrix(0, 1, d)
  beta_SE <- matrix(0, 1, d)
  P_raw_DLASSO <- matrix(0, 1, d)
  M1 <- t(t(M_raw[, 1]))

  checkParallel("hima_microbiome", parallel, ncore, verbose)
  
  if (verbose) message("Step 1: ILR Transformation and De-biased Lasso estimates ...", "  (", format(Sys.time(), "%X"), ")")

  if (verbose) {
    if (is.null(COV)) {
      message("        No covariate was adjusted.")
    } else {
      message("        Adjusting for covariate(s): ", paste0(colnames(COV), collapse = ", "))
    }
  }
  
  results_loop <- foreach(k = seq_len(d), .combine = rbind) %dopar% {
    M <- M_raw
    M[, 1] <- M[, k]
    M[, k] <- M1
    MT <- matrix(0, n, d - 1)
    for (i in 1:n) {
      for (j in 1:(d - 1)) {
        C_1 <- sqrt((d - j) / (d - j + 1))
        C_2 <- prod(M[i, (j + 1):d]^(1 / (d - j)))
        MT[i, j] <- C_1 * log(M[i, j] / C_2)
      }
    }

    MT <- scale(MT)
    MX <- cbind(MT, X)

    fit.dlasso <- DLASSO_fun(MX, Y)

    beta_est <- fit.dlasso[1]
    beta_se <- fit.dlasso[2]
    P_b <- 2 * (1 - pnorm(abs(beta_est / beta_se), 0, 1))

    lm.fit <- stats::lm(MT[, 1] ~ X)
    lm.out <- summary(lm.fit)
    alpha_est <- lm.out$coefficients[2, 1]
    alpha_se <- lm.out$coefficients[2, 2]
    P_a <- 2 * (1 - pnorm(abs(alpha_est / alpha_se), 0, 1))
    c(beta_est, beta_se, alpha_est, alpha_se, max(P_a, P_b))
  }
  if (is.null(dim(results_loop))) results_loop <- matrix(results_loop, nrow = 1)
  beta_EST <- results_loop[, 1]
  beta_SE <- results_loop[, 2]
  alpha_EST <- results_loop[, 3]
  alpha_SE <- results_loop[, 4]
  P_raw_DLASSO <- results_loop[, 5]
  # the end of k
  
  P_adj_DLASSO <- as.numeric(P_raw_DLASSO)
  
  if (verbose) message("Step 2: Closed testing-based procedure ...", "     (", format(Sys.time(), "%X"), ")")
  
  ## The FDR method
  set <- which(P_adj_DLASSO < FDRcut)
  hom <- hommel::hommel(P_adj_DLASSO, simes = FALSE)
  N1 <- hommel::discoveries(hom, set, incremental = TRUE, alpha = FDRcut)
  
  if (length(set) > 0) {
    L <- length(set)
    N2 <- matrix(0, 1, L)
    N2[2:L] <- N1[1:(L - 1)]
  }
  
  N0 <- N1 - N2
  
  ID_FDR <- set[which(N0 > 0)]
  
  IDE <- alpha_EST[ID_FDR] * beta_EST[ID_FDR]
  
  if (length(ID_FDR) > 0) {
    out_result <- data.frame(
      Index = M_ID_name[ID_FDR],
      alpha_hat = alpha_EST[ID_FDR],
      alpha_se = alpha_SE[ID_FDR],
      beta_hat = beta_EST[ID_FDR],
      beta_se = beta_SE[ID_FDR],
      IDE = IDE,
      rimp = abs(IDE) / sum(abs(IDE)) * 100,
      pmax = P_adj_DLASSO[ID_FDR], row.names = NULL
    )
    if (verbose) message(paste0("        ", length(ID_FDR), " significant mediator(s) identified."))
  } else {
    if (verbose) message("        No significant mediator identified.")
    out_result <- NULL
  }
  
  if (verbose) message("Done!", "     (", format(Sys.time(), "%X"), ")")

  doParallel::stopImplicitCluster()

  return(out_result)
}
