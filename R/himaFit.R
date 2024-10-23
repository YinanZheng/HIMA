# This is the wrapper function for high-dimensional mediation analysis
#' High-dimensional Mediation Analysis
#'
#' \code{himaFit} is a wrapper function designed to perform various HIMA methods for estimating and testing high-dimensional mediation effects.
#' \code{himaFit} can automatically select the appropriate HIMA method based on the outcome and mediator data type.
#'
#' @param formula an object of class \code{formula}: a symbolic description of the overall effect model, i.e.,
#' \code{outcome ~ exposure + covariates}, to be fitted. Make sure the "exposure" is the variable of interest, which
#' must be listed as the first variable in the right hand side of the formula.
#' @param data.pheno a data frame containing exposure and covariates that are listed in the right hand side of the \code{formula}.
#' The variable names must match those listed in \code{formula}. By default \code{himaFit} will scale \code{data.pheno}.
#' @param data.M a \code{data.frame} or \code{matrix} of high-dimensional mediators. Rows represent samples, columns
#' represent variables. By default \code{himaFit} will scale \code{data.M}.
#' @param mediator.type data type of high-dimensional mediators (\code{data.M}). Either \code{'gaussian'} (default, for continuous mediators),
#' \code{'negbin'} (i.e., negative binomial, for RNA-seq data as mediators), or \code{'compositional'} (for microbiome data as mediators).
#' @param penalty the penalty to be applied to the model. Either \code{'DBlasso'} (De-biased LASSO, default),
#' \code{'MCP'}, \code{'SCAD'}, or \code{'lasso'}. Please note, survival HIMA and microbiome HIMA can be only performed with \code{'DBlasso'};
#' Quantile HIMA cannot be performed with \code{'DBlasso'}. Not applicable for efficient HIMA (when \code{efficient = TRUE}), as it will always apply \code{'MCP'}.
#' @param quantile use quantile HIMA (\code{qHIMA}). Only applicable for classic HIMA with continues outcome. Please add parameter \code{tau} to specify the
#' desired quantile (can be a vector), otherwise the default \code{tau = 0.5} will be used.
#' @param efficient use efficient HIMA (\code{eHIMA}). Only applicable for classic HIMA with continuous outcome.
#' and \code{mediator.type = "gaussian"}). Default = \code{FALSE}.
#' @param scale logical. Should the function scale the data (exposure, mediators, and covariates)? Default = \code{TRUE}. Note: for simulation study, scale
#' can be turned off to avoid estimate compression.
#' @param Sigcut cutoff applied to select significant mediators. Default = \code{0.05}.
#' @param verbose logical. Should the function be verbose and show the progression? Default = \code{FALSE}.
#' @param ... reserved passing parameter.
#'
#' @return A data.frame containing mediation testing results of selected mediators.
#' \describe{
#'     \item{ID: }{Mediator ID/name.}
#'     \item{alpha: }{Coefficient estimates of exposure (X) --> mediators (M) (adjusted for covariates).}
#'     \item{beta: }{Coefficient estimates of mediators (M) --> outcome (Y) (adjusted for covariates and exposure).}
#'     \item{alpha*beta: }{Mediation (indirect) effect.}
#'     \item{Relative Importance: }{Relative importance of the mediator. It is the proportion of the mediation effect for each mediator
#'     out of the sum of the mediation effect (absolute value) across all significant mediators selected.}
#'     \item{p-value: }{Joint raw p-value of significant mediators selected based on corresponding approach.}
#'     \item{tau: }{Quantile level of the outcome (applicable only to the quantile mediation model).}
#' }
#'
#' @references
#' 1. Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J, Just A, Colicino E, Vokonas P, Zhao L,
#' Lv J, Baccarelli A, Hou L, Liu L. Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies.
#' Bioinformatics. 2016. DOI: 10.1093/bioinformatics/btw351. PMID: 27357171; PMCID: PMC5048064
#'
#' 2. Zhang H, Zheng Y, Hou L, Zheng C, Liu L. Mediation Analysis for Survival Data with High-Dimensional Mediators.
#' Bioinformatics. 2021. DOI: 10.1093/bioinformatics/btab564. PMID: 34343267; PMCID: PMC8570823
#'
#' 3. Zhang H, Chen J, Feng Y, Wang C, Li H, Liu L. Mediation Effect Selection in High-dimensional and Compositional Microbiome data.
#' Stat Med. 2021. DOI: 10.1002/sim.8808. PMID: 33205470; PMCID: PMC7855955
#'
#' 4. Zhang H, Chen J, Li Z, Liu L. Testing for Mediation Effect with Application to Human Microbiome Data.
#' Stat Biosci. 2021. DOI: 10.1007/s12561-019-09253-3. PMID: 34093887; PMCID: PMC8177450
#'
#' 5. Perera C, Zhang H, Zheng Y, Hou L, Qu A, Zheng C, Xie K, Liu L. HIMA2: High-dimensional Mediation Analysis and Its Application in
#' Epigenome-wide DNA Methylation Data. BMC Bioinformatics. 2022. DOI: 10.1186/s12859-022-04748-1. PMID: 35879655; PMCID: PMC9310002
#'
#' 6. Zhang H, Hong X, Zheng Y, Hou L, Zheng C, Wang X, Liu L. High-Dimensional Quantile Mediation Analysis with Application to a Birth
#' Cohort Study of Mother–Newborn Pairs. Bioinformatics. 2024. DOI: 10.1093/bioinformatics/btae055. PMID: 38290773; PMCID: PMC10873903
#'
#' 7. Bai X, Zheng Y, Hou L, Zheng C, Liu L, Zhang H. An Efficient Testing Procedure for High-dimensional Mediators with FDR Control.
#' Statistics in Biosciences. 2024. DOI: 10.1007/s12561-024-09447-4.
#'
#' @examples
#' \dontrun{
#' # Note: In the following examples, M1, M2, and M3 are true mediators.
#' data(himaDat)
#'
#' # Example 1 (continuous outcome - linear HIMA):
#' head(himaDat$Example1$PhenoData)
#'
#' e1 <- himaFit(Outcome ~ Treatment + Sex + Age,
#'   data.pheno = himaDat$Example1$PhenoData,
#'   data.M = himaDat$Example1$Mediator,
#'   mediator.type = "gaussian",
#'   penalty = "MCP", # Can be "DBlasso" for dblassoHIMA
#'   scale = FALSE
#' ) # Disabled only for simulation data
#' e1
#' attributes(e1)$variable.labels
#'
#' # Efficient HIMA (only applicable to mediators and outcomes that are
#' # both continuous and normally distributed.)
#' e1e <- himaFit(Outcome ~ Treatment + Sex + Age,
#'   data.pheno = himaDat$Example1$PhenoData,
#'   data.M = himaDat$Example1$Mediator,
#'   mediator.type = "gaussian",
#'   efficient = TRUE,
#'   penalty = "MCP", # Efficient HIMA does not support DBlasso
#'   scale = FALSE
#' ) # Disabled only for simulation data
#' e1e
#' attributes(e1e)$variable.labels
#'
#' # Example 2 (binary outcome - logistic HIMA):
#' head(himaDat$Example2$PhenoData)
#'
#' e2 <- himaFit(Disease ~ Treatment + Sex + Age,
#'   data.pheno = himaDat$Example2$PhenoData,
#'   data.M = himaDat$Example2$Mediator,
#'   mediator.type = "gaussian",
#'   penalty = "MCP",
#'   scale = FALSE
#' ) # Disabled only for simulation data
#' e2
#' attributes(e2)$variable.labels
#'
#' # Example 3 (time-to-event outcome - survival HIMA):
#' head(himaDat$Example3$PhenoData)
#'
#' e3 <- himaFit(Surv(Status, Time) ~ Treatment + Sex + Age,
#'   data.pheno = himaDat$Example3$PhenoData,
#'   data.M = himaDat$Example3$Mediator,
#'   mediator.type = "gaussian",
#'   penalty = "DBlasso",
#'   scale = FALSE
#' ) # Disabled only for simulation data
#' e3
#' attributes(e3)$variable.labels
#'
#' # Example 4 (compositional data as mediator, e.g., microbiome):
#' head(himaDat$Example4$PhenoData)
#'
#' e4 <- himaFit(Outcome ~ Treatment + Sex + Age,
#'   data.pheno = himaDat$Example4$PhenoData,
#'   data.M = himaDat$Example4$Mediator,
#'   mediator.type = "compositional",
#'   penalty = "DBlasso"
#' ) # Scaling is always enabled internally for microHIMA
#' e4
#' attributes(e4)$variable.labels
#'
#' #' # Example 5 (quantile mediation anlaysis - quantile HIMA):
#' head(himaDat$Example5$PhenoData)
#'
#' # Note that the function will prompt input for quantile level.
#' e5 <- himaFit(Outcome ~ Treatment + Sex + Age,
#'   data.pheno = himaDat$Example5$PhenoData,
#'   data.M = himaDat$Example5$Mediator,
#'   mediator.type = "gaussian",
#'   quantile = TRUE,
#'   penalty = "MCP", # Quantile HIMA does not support DBlasso
#'   scale = FALSE, # Disabled only for simulation data
#'   tau = c(0.3, 0.5, 0.7)
#' ) # Specify multiple quantile level
#' e5
#' attributes(e5)$variable.labels
#' }
#'
#' @export
himaFit <- function(formula,
                 data.pheno,
                 data.M,
                 mediator.type = c("gaussian", "negbin", "compositional"),
                 penalty = c("DBlasso", "MCP", "SCAD", "lasso"),
                 quantile = FALSE,
                 efficient = FALSE,
                 scale = TRUE,
                 Sigcut = 0.05,
                 verbose = FALSE,
                 ...) {
  # extract data from formula and input dataset
  d <- .f_extract(data.pheno = data.pheno, f = formula)

  if (d$type == "categorical") stop("Currently HIMA does not support categorical outcome.")

  if (efficient && quantile) stop("Cannot have both 'efficient' and 'quantile' set to TRUE simultaneously.")

  mediator.type <- match.arg(mediator.type)
  penalty <- match.arg(penalty)

  if (Sigcut > 1) Sigcut <- 1

  if (is.null(colnames(data.M))) colnames(data.M) <- seq_len(ncol(data.M))

  # Conditions where 'DBlasso' is not applicable
  invalid_dblasso <- efficient || quantile || d$type == "binary" || mediator.type == "negbin"

  # Conditions where 'DBlasso' must be used
  require_dblasso <- d$type == "survival" || mediator.type == "compositional"

  # Check for conflicting conditions
  if (invalid_dblasso && require_dblasso) {
    stop("Conflicting conditions: Cannot determine appropriate penalty under the current settings.")
  }

  # Adjust 'penalty' based on the conditions
  if (penalty == "DBlasso" && invalid_dblasso) {
    message("Note: The selected conditions do not support de-biased Lasso penalty. Switching to 'MCP' ...")
    penalty <- "MCP"
  } else if (penalty != "DBlasso" && require_dblasso) {
    message("Note: The selected conditions require de-biased Lasso penalty. Switching to 'DBlasso' ...")
    penalty <- "DBlasso"
  }

  ######

  # eHIMA / qHIMA
  if (efficient || quantile) {
    if (efficient) {
      message("Running efficient HIMA with ", penalty, " penalty...")
      if (d$type != "continuous" || mediator.type != "gaussian") {
        stop("Efficient HIMA is only applicable to mediators and outcomes that are BOTH continuous and normally distributed.")
      }

      res <- eHIMA(
        X = d$X,
        M = data.M,
        Y = d$Y,
        COV = d$COV,
        topN = NULL,
        scale = scale,
        FDRcut = Sigcut,
        verbose = verbose
      )

      results <- .res_prep(res, method_text = "DACT method with BH-FDR", Sigcut = Sigcut)
    }

    if (quantile) {
      message("Running quantile HIMA with ", penalty, " penalty...")
      if (d$type != "continuous" || mediator.type != "gaussian") {
        stop("Quantile HIMA is only applicable to mediators and outcomes that are BOTH continuous and normally distributed.")
      }

      # tau <- readline(prompt = "Enter quantile level(s) (between 0-1, multiple values accepted): ")
      # tau <- eval(parse(text = paste0("c(", tau, ")")))
      res <- qHIMA(
        X = d$X,
        M = data.M,
        Y = d$Y,
        COV = d$COV,
        penalty = penalty,
        topN = NULL,
        scale = scale,
        Bonfcut = 0.05,
        verbose = verbose,
        ...
      )

      results <- .res_prep(res, method_text = "Bonferroni-adjusted p", Sigcut = Sigcut, q = TRUE)
    }
  } else { # If not eHIMA AND not qHIMA
    # DBlasso
    if (penalty == "DBlasso") {
      if (d$type == "continuous") {
        if (mediator.type == "gaussian") {
          message("Running linear HIMA with de-biased Lasso penalty...")
          res <- dblassoHIMA(
            X = d$X,
            M = data.M,
            Y = d$Y,
            COV = d$COV,
            topN = NULL,
            scale = scale,
            FDRcut = Sigcut,
            verbose = verbose
          )

          results <- .res_prep(res, method_text = "HDMT pointwise FDR", Sigcut = Sigcut)
        } else if (mediator.type == "compositional") {
          message("Running compositional HIMA with de-biased Lasso penalty...")
          res <- microHIMA(
            X = d$X,
            OTU = data.M,
            Y = d$Y,
            COV = d$COV,
            FDRcut = Sigcut,
            verbose = verbose
          )

          results <- .res_prep(res, method_text = "Hommel FDR", Sigcut = Sigcut)
        }
      } else if (d$type == "survival") {
        message("Running survival HIMA with de-biased Lasso penalty...")
        res <- survHIMA(
          X = d$X,
          M = data.M,
          OT = d$Y$OT,
          status = d$Y$status,
          COV = d$COV,
          topN = NULL,
          scale = scale,
          FDRcut = Sigcut,
          verbose = verbose
        )

        results <- .res_prep(res, method_text = "HDMT pointwise FDR", Sigcut = Sigcut)
      }
    } else { # If penalty is not DBlasso
      res <- classicHIMA(
        X = d$X,
        M = data.M,
        Y = d$Y,
        COV.XM = d$COV,
        Y.type = d$type,
        M.type = mediator.type,
        penalty = penalty,
        topN = NULL,
        parallel = FALSE,
        ncore = 1,
        scale = scale,
        Bonfcut = Sigcut,
        verbose = verbose
      )

      results <- .res_prep(res, method_text = "Bonferroni-adjusted p", Sigcut = Sigcut)
    }
  }

  if (is.null(results)) message("No significant mediator found!")

  return(results)
}


### Internal function: Result prepare function
.res_prep <- function(res, method_text, Sigcut, q = FALSE) {
  if (is.null(res)) {
    return(NULL)
  }
  results <- data.frame(
    ID = res$Index,
    alpha = res$alpha_hat,
    beta = res$beta_hat,
    `alpha*beta` = res$IDE,
    `Relative Importance (%)` = res$rimp,
    `p-value` = res$pmax,
    check.names = FALSE
  )
  results <- results[order(results$`p-value`), ]
  attr_test <- c(
    "ID: Mediator ID",
    "alpha: Effect of exposure on mediator",
    "beta: Effect of mediator on outcome",
    "alpha*beta: Mediation (indirect) effect",
    "Relative Importance (%): Relative importance of the mediator out of all significant mediators",
    paste0("p-value: Joint raw p-value of significant mediator selected based on ", method_text, " < ", Sigcut)
  )
  if (q) {
    results <- data.frame(results,
      tau = res$tau,
      check.names = FALSE
    )
    attr_test <- c(
      attr_test,
      "tau: Quantile level of the outcome"
    )
  }
  attr(results, "variable.labels") <- attr_test
  rownames(results) <- NULL
  return(results)
}



### Internal function: recognize character variables and convert to dummy (only in himaFit)
.convert_to_dummies <- function(df) {
  char_cols <- sapply(df, is.character)
  df[char_cols] <- lapply(df[char_cols], as.factor)
  df_dummies <- as.data.frame(model.matrix(~ . - 1, data = df))
  non_char_cols <- names(df)[!char_cols]
  df_dummies[non_char_cols] <- df[non_char_cols]
  return(df_dummies)
}



### Internal function: determine_data_type: A function to determine data type
.data_type <- function(x, categorical_threshold = 10) {
  # Remove NA values
  x_clean <- x[!is.na(x)]
  # Get unique values
  uniq_vals <- unique(x_clean)
  num_unique <- length(uniq_vals)

  # Handle variables with no data
  if (num_unique == 0) {
    stop("No valid outcome data detected!") # Variable contains only NA values
  }

  # Classify binary variables
  if (num_unique == 1 || num_unique == 2) {
    return("binary")
  }

  # Classify continuous and categorical variables
  if (is.numeric(x) || is.integer(x)) {
    if (num_unique > categorical_threshold) {
      return("continuous")
    } else {
      return("categorical")
    }
  } else if (is.factor(x) || is.character(x)) {
    return("categorical")
  } else if (is.logical(x)) {
    return("binary")
  } else {
    return("categorical") # Default classification
  }
}



### Internal function: Formula extractor function
.f_extract <- function(data.pheno, f) {
  if (as.character(f[[2]])[1] == "Surv") {
    response_vars <- as.character(f[[2]])[c(2, 3)]
    ind_vars <- all.vars(f)[-c(1, 2)]
    Y <- list(
      status = data.pheno[, response_vars[1]],
      OT = data.pheno[, response_vars[2]]
    )
    type <- "survival"
  } else {
    response_var <- as.character(f[[2]])
    ind_vars <- all.vars(f)[-1]
    Y <- data.pheno[, response_var]
    type <- .data_type(Y)
  }

  X <- data.pheno[, ind_vars[1]]

  if (length(ind_vars) > 1) {
    COV <- .convert_to_dummies(data.pheno[, ind_vars[-1]])
  } else {
    COV <- NULL
  }

  return(list(
    Y = Y,
    X = X,
    COV = COV,
    type = type
  ))
}



### hima2 function has been renamed to himaFit with parameters updated.
hima2 <- function(formula,
                  data.pheno,
                  data.M,
                  outcome.family = c("gaussian", "binomial", "survival", "quantile"),
                  mediator.family = c("gaussian", "negbin", "compositional"),
                  penalty = c("DBlasso", "MCP", "SCAD", "lasso"),
                  efficient = FALSE,
                  scale = TRUE,
                  Sigcut = 0.05,
                  verbose = FALSE,
                  ...) {
  .Deprecated(msg = "hima2() is deprecated. Please use himaFit() with updated parameters.")

  outcome.family <- match.arg(outcome.family)
  if (outcome.family == "quantile") quantile <- TRUE else quantile <- FALSE
  if (!missing(outcome.family)) {
    warning("'outcome.family' is deprecated. The outcome data type can now be automatically recognized.")
  }

  mediator.type <- match.arg(mediator.family)
  penalty <- match.arg(penalty)

  himaFit(formula,
    data.pheno,
    data.M,
    mediator.type = mediator.type,
    penalty = penalty,
    quantile = quantile,
    efficient = efficient,
    scale = scale,
    Sigcut = Sigcut,
    verbose = verbose,
    ...
  )
}
