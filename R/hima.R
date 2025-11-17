# This is the wrapper function for high-dimensional mediation analysis
#' High-dimensional Mediation Analysis
#'
#' \code{hima} is a wrapper function designed to perform various HIMA methods for estimating and testing high-dimensional mediation effects.
#' \code{hima} can automatically select the appropriate HIMA method based on the outcome and mediator data type.
#'
#' @param formula an object of class \code{formula} representing the overall effect model to be fitted, specified as \code{outcome ~ exposure + covariates}.
#' The "exposure" variable (the variable of interest) must be listed first on the right-hand side of the formula.
#' For survival outcomes specified using \code{Surv()}, the exposure should be the first
#' variable after the \code{\link[base:tilde]{~}}.
#' @param data.pheno a data frame containing the exposure, outcome, and covariates specified in the formula. Variable names in \code{data.pheno} must match those 
#' in the formula. When \code{scale = TRUE}, the exposure and covariates will be scaled (the outcome retains its original scale).
#' @param data.M a \code{data.frame} or \code{matrix} of high-dimensional mediators, with rows representing samples and columns representing mediator variables. 
#' When \code{scale = TRUE}, \code{data.M} will be scaled.
#' @param mediator.type a character string indicating the data type of the high-dimensional mediators (\code{data.M}). Options are: \code{'gaussian'} (default): 
#' for continuous mediators. \code{'negbin'}: for count data mediators modeled using the negative binomial distribution (e.g., RNA-seq data). \code{'compositional'}: 
#' for compositional data mediators (e.g., microbiome data).
#' @param penalty a character string specifying the penalty method to apply in the model. Options are: \code{'DBlasso'}: De-biased LASSO (default). \code{'MCP'}: 
#' Minimax Concave Penalty. \code{'SCAD'}: Smoothly Clipped Absolute Deviation. \code{'lasso'}: Least Absolute Shrinkage and Selection Operator. Note: Survival HIMA and microbiome 
#' HIMA can only be performed with \code{'DBlasso'}. Quantile HIMA and efficient HIMA cannot use \code{'DBlasso'}; they always apply \code{'MCP'}.
#' @param quantile logical. Indicates whether to use quantile HIMA (\code{hima_quantile}). Default is \code{FALSE}. Applicable only for classic HIMA with a continuous outcome and 
#' \code{mediator.type = 'gaussian'}. If \code{TRUE}, specify the desired quantile(s) using the \code{tau} parameter; otherwise, the default \code{tau = 0.5} (i.e., median) is used.
#' @param efficient logical. Indicates whether to use efficient HIMA (\code{hima_efficient}). Default is \code{FALSE}. Applicable only for classic HIMA with a continuous outcome and 
#' \code{mediator.type = 'gaussian'}.
#' @param longitudinal logical. Indicates whether to run the longitudinal survival mediation model \code{hima_survival_long} (requires a \code{Surv(tstart, tstop, status)} outcome). 
#' Default = \code{FALSE}.
#' @param id.var Character string specifying the column name in \code{data.pheno} that stores subject identifiers. Required when \code{longitudinal = TRUE}.
#' @param scale logical. Determines whether the function scales the data (exposure, mediators, and covariates). Default is \code{TRUE}. Note: For simulation studies, set 
#' \code{scale = FALSE} to avoid estimate compression (i.e., shrinkage of estimates toward zero due to scaling).
#' @param sigcut numeric. The significance cutoff for selecting mediators. Default is \code{0.05}.
#' @param contrast a named list of contrasts to be applied to factor variables in the covariates (cannot be the variable of interest).
#' @param subset an optional vector specifying a subset of observations to use in the analysis.
#' @param verbose logical. Determines whether the function displays progress messages. Default is \code{FALSE}.
#' @param parallel logical. Enable parallel computing feature? Default = \code{FALSE}.
#' @param ncore number of cores to run parallel computing Valid when \code{parallel = TRUE}.
#' @param ... reserved passing parameter (or for future use).
#'
#' @return A data.frame containing mediation testing results of selected mediators.
#' \describe{
#'     \item{ID: }{Mediator ID/name.}
#'     \item{alpha: }{Coefficient estimates of exposure (X) --> mediators (M) (adjusted for covariates).}
#'     \item{beta: }{Coefficient estimates of mediators (M) --> outcome (Y) (adjusted for covariates and exposure).}
#'     \item{alpha*beta: }{The estimated indirect (mediation) effect of exposure on outcome through each mediator.}
#'     \item{rimp: }{Relative importance- the proportion of each mediator's mediation effect relative to the sum of the absolute mediation effects of all significant mediators.}
#'     \item{p-value: }{The joint p-value assessing the significance of each mediator's indirect effect, calculated based on the corresponding statistical approach.}
#'     \item{tau: }{The quantile level of the outcome (applicable only when using the quantile mediation model).}
#' }
#'
#' @references
#' 1. Zhang H, Zheng Y, Hou L, Liu L, HIMA: An R Package for High-Dimensional Mediation Analysis. Journal of Data Science. 2025. DOI: 10.6339/25-JDS1192
#' 
#' 2. Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J, Just A, Colicino E, Vokonas P, Zhao L,
#' Lv J, Baccarelli A, Hou L, Liu L. Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies.
#' Bioinformatics. 2016. DOI: 10.1093/bioinformatics/btw351. PMID: 27357171; PMCID: PMC5048064
#' 
#' 3. Zhang H, Zheng Y, Hou L, Zheng C, Liu L. Mediation Analysis for Survival Data with High-Dimensional Mediators.
#' Bioinformatics. 2021. DOI: 10.1093/bioinformatics/btab564. PMID: 34343267; PMCID: PMC8570823
#' 
#' 4. Zhang H, Chen J, Feng Y, Wang C, Li H, Liu L. Mediation Effect Selection in High-dimensional and Compositional Microbiome data.
#' Stat Med. 2021. DOI: 10.1002/sim.8808. PMID: 33205470; PMCID: PMC7855955
#' 
#' 5. Zhang H, Chen J, Li Z, Liu L. Testing for Mediation Effect with Application to Human Microbiome Data.
#' Stat Biosci. 2021. DOI: 10.1007/s12561-019-09253-3. PMID: 34093887; PMCID: PMC8177450
#' 
#' 6. Perera C, Zhang H, Zheng Y, Hou L, Qu A, Zheng C, Xie K, Liu L. HIMA2: High-dimensional Mediation Analysis and Its Application in
#' Epigenome-wide DNA Methylation Data. BMC Bioinformatics. 2022. DOI: 10.1186/s12859-022-04748-1. PMID: 35879655; PMCID: PMC9310002
#' 
#' 7. Zhang H, Hong X, Zheng Y, Hou L, Zheng C, Wang X, Liu L. High-Dimensional Quantile Mediation Analysis with Application to a Birth
#' Cohort Study of Mother–Newborn Pairs. Bioinformatics. 2024. DOI: 10.1093/bioinformatics/btae055. PMID: 38290773; PMCID: PMC10873903
#' 
#' 8. Bai X, Zheng Y, Hou L, Zheng C, Liu L, Zhang H. An Efficient Testing Procedure for High-dimensional Mediators with FDR Control.
#' Statistics in Biosciences. 2024. DOI: 10.1007/s12561-024-09447-4.
#' 
#' 9. Liu L, Zhang H, Zheng Y, Gao T, Zheng C, Zhang K, Hou L, Liu L. High-dimensional mediation analysis for longitudinal mediators and 
#' survival outcomes. Briefings in Bioinformatics. 2025. DOI: 10.1093/bib/bbaf206. PMID: 40350699 PMCID: PMC12066418
#'
#' @examples
#' \dontrun{
#' # Note: In the following examples, M1, M2, and M3 are true mediators.
#'
#' # Example 1 (continuous outcome - linear HIMA):
#' data(ContinuousOutcome)
#' pheno_data <- ContinuousOutcome$PhenoData
#' mediator_data <- ContinuousOutcome$Mediator
#'
#' e1 <- hima(Outcome ~ Treatment + Sex + Age,
#'   data.pheno = pheno_data,
#'   data.M = mediator_data,
#'   mediator.type = "gaussian",
#'   penalty = "MCP", # Can be "DBlasso" for hima_dblasso
#'   scale = FALSE, # Disabled only for simulation data
#'   verbose = TRUE
#' ) 
#' summary(e1)
#'
#' # Efficient HIMA (only applicable to mediators and outcomes that are
#' # both continuous and normally distributed.)
#' e1e <- hima(Outcome ~ Treatment + Sex + Age,
#'   data.pheno = pheno_data,
#'   data.M = mediator_data,
#'   mediator.type = "gaussian",
#'   efficient = TRUE,
#'   penalty = "MCP", # Efficient HIMA does not support DBlasso
#'   scale = FALSE, # Disabled only for simulation data
#'   verbose = TRUE
#' ) 
#' summary(e1e)
#'
#' # Example 2 (binary outcome - logistic HIMA):
#' data(BinaryOutcome)
#' pheno_data <- BinaryOutcome$PhenoData
#' mediator_data <- BinaryOutcome$Mediator
#'
#' e2 <- hima(Disease ~ Treatment + Sex + Age,
#'   data.pheno = pheno_data,
#'   data.M = mediator_data,
#'   mediator.type = "gaussian",
#'   penalty = "MCP",
#'   scale = FALSE, # Disabled only for simulation data
#'   verbose = TRUE
#' ) 
#' summary(e2)
#'
#' # Example 3 (time-to-event outcome - survival HIMA):
#' data(SurvivalData)
#' pheno_data <- SurvivalData$PhenoData
#' mediator_data <- SurvivalData$Mediator
#'
#' e3 <- hima(Surv(Time, Status) ~ Treatment + Sex + Age,
#'   data.pheno = pheno_data,
#'   data.M = mediator_data,
#'   mediator.type = "gaussian",
#'   penalty = "DBlasso",
#'   scale = FALSE, # Disabled only for simulation data
#'   verbose = TRUE
#' ) # Parallel computing feature is recommended
#' summary(e3)
#' 
#' # Longitudinal mediator + survival HIMA:
#' data(SurvivalLongData)
#' pheno_data <- SurvivalLongData$PhenoData
#' mediator_data <- SurvivalLongData$Mediator
#'
#' e3long <- hima(Surv(Tstart, Tstop, Status) ~ Treatment + Sex + Age,
#'   data.pheno = pheno_data,
#'   data.M = mediator_data,
#'   mediator.type = "gaussian",
#'   penalty = "lasso",
#'   longitudinal = TRUE,
#'   id.var = "ID",
#'   scale = FALSE, # Disabled only for simulation data
#'   verbose = TRUE
#' ) # Parallel computing feature is recommended
#' summary(e3long)
#' 
#' # Example 4 (compositional data as mediator, e.g., microbiome):
#' data(MicrobiomeData)
#' pheno_data <- MicrobiomeData$PhenoData
#' mediator_data <- MicrobiomeData$Mediator
#'
#' e4 <- hima(Outcome ~ Treatment + Sex + Age,
#'   data.pheno = pheno_data,
#'   data.M = mediator_data,
#'   mediator.type = "compositional",
#'   penalty = "DBlasso",
#'   verbose = TRUE
#' ) # Scaling is always enabled internally for hima_microbiome
#' summary(e4)
#'
#' #' # Example 5 (quantile mediation analysis - quantile HIMA):
#' data(QuantileData)
#' pheno_data <- QuantileData$PhenoData
#' mediator_data <- QuantileData$Mediator
#'
#' # Note that the function will prompt input for quantile level.
#' e5 <- hima(Outcome ~ Treatment + Sex + Age,
#'   data.pheno = pheno_data,
#'   data.M = mediator_data,
#'   mediator.type = "gaussian",
#'   quantile = TRUE,
#'   penalty = "MCP", # Quantile HIMA does not support DBlasso
#'   scale = FALSE, # Disabled only for simulation data
#'   tau = c(0.3, 0.5, 0.7),
#'   verbose = TRUE
#' ) # Specify multiple quantile level
#' summary(e5)
#' }
#'
#' @export
hima <- function(formula,
                 data.pheno,
                 data.M,
                 mediator.type = c("gaussian", "negbin", "compositional"),
                 penalty = c("DBlasso", "MCP", "SCAD", "lasso"),
                 quantile = FALSE,
                 efficient = FALSE,
                 longitudinal = FALSE,
                 id.var = NULL,
                 scale = TRUE,
                 sigcut = 0.05,
                 contrast = NULL,
                 subset = NULL,
                 verbose = FALSE,
                 parallel = FALSE,
                 ncore = 1,
                 ...) {
  
  # Handle subset
  if (!is.null(subset)) {
    subset_rows <- eval(subset, data.pheno, parent.frame())
    if (!is.logical(subset_rows) || length(subset_rows) != nrow(data.pheno)) {
      stop("Subset must be a logical vector with length equal to the number of rows \nin data.pheno.")
    }
    data.pheno <- data.pheno[subset_rows, , drop = FALSE]
    data.M <- data.M[subset_rows, , drop = FALSE]
  }
  
  # Identify the variable of interest (first term on the right-hand side of the formula)
  if (as.character(formula[[2]])[1] == "Surv") {
    rhs_vars <- all.vars(formula[[3]])
    var_of_interest <- rhs_vars[1]
  } else {
    var_of_interest <- all.vars(formula)[2]
  }
  
  # Apply contrasts only to covariates (not to the variable of interest)
  if (!is.null(contrast)) {
    for (var in names(contrast)) {
      if (var == var_of_interest) {
        warning(sprintf("Contrast is ignored for the variable of interest: '%s'.", var))
        next
      }
      if (var %in% names(data.pheno)) {
        if (is.factor(data.pheno[[var]])) {
          contrasts(data.pheno[[var]]) <- contrast[[var]]
        } else {
          warning(sprintf("Contrast not applied: Variable '%s' is not a factor.", var))
        }
      } else {
        warning(sprintf("Variable '%s' specified in contrasts does not exist in data.pheno.", var))
      }
    }
  }
  
  # extract data from formula and input dataset
  d <- .f_extract(data.pheno = data.pheno, f = formula)
  
  # check mutually exclusive analysis modes
  mode_flags <- c(
    efficient   = isTRUE(efficient),
    quantile    = isTRUE(quantile),
    longitudinal = isTRUE(longitudinal)
  )
  if (sum(mode_flags) > 1L) {
    stop("Arguments 'efficient', 'quantile', and 'longitudinal' are mutually exclusive; please set at most one of them to TRUE.")
  }
  
  if (d$type == "categorical") stop("Currently HIMA does not support categorical outcome.")
  
  if (longitudinal && d$type != "survival") {
    stop("Longitudinal mediation is only available for survival outcomes specified via Surv().")
  }
  
  id_vec <- NULL
  
  if (longitudinal) {
    if (is.null(d$Y$tstart) || is.null(d$Y$tstop)) {
      stop("Please specify the survival outcome as Surv(tstart, tstop, status) when longitudinal = TRUE.")
    }
    if (is.null(id.var)) {
      stop("Argument 'id.var' must be provided when longitudinal = TRUE.")
    }
    if (!id.var %in% names(data.pheno)) {
      stop("The column specified in 'id.var' does not exist in data.pheno.")
    }
    id_vec <- data.pheno[[id.var]]
  }
  
  mediator.type <- match.arg(mediator.type)
  penalty <- match.arg(penalty)
  
  # Longitudinal: penalty must be lasso
  if (longitudinal && penalty != "lasso") {
    if (verbose) message("Note: For longitudinal survival mediation, only 'lasso' penalty is supported. \nSwitching to 'lasso' ...")
    penalty <- "lasso"
  }
  
  if (sigcut > 1) sigcut <- 1
  
  if (is.null(colnames(data.M))) colnames(data.M) <- seq_len(ncol(data.M))
  
  # Conditions where 'DBlasso' is not applicable
  invalid_dblasso <- efficient || quantile || longitudinal || d$type == "binary" || mediator.type == "negbin"
  
  # Conditions where 'DBlasso' must be used
  require_dblasso <- (!longitudinal) && (d$type == "survival" || mediator.type == "compositional")
  
  # Check for conflicting conditions
  if (invalid_dblasso && require_dblasso) {
    stop("Conflicting conditions: Cannot determine appropriate penalty under the current settings.")
  }
  
  # Adjust 'penalty' based on the conditions
  if (penalty == "DBlasso" && invalid_dblasso) {
    if (verbose) message("Note: The selected conditions do not support de-biased Lasso penalty. \nSwitching to 'MCP' ...")
    penalty <- "MCP"
  } else if (penalty != "DBlasso" && require_dblasso) {
    if (verbose) message("Note: The selected conditions require de-biased Lasso penalty. \nSwitching to 'DBlasso' ...")
    penalty <- "DBlasso"
  }
  
  ######
  
  # hima_efficient / hima_quantile / hima_survival_long
  if (efficient || quantile || longitudinal) {
    if (efficient) {
      if (penalty != "MCP") {
        if (verbose) message("Note: Efficient HIMA works best with 'MCP' penalty. Switching penalty to 'MCP'...")
      }
      if (verbose) message("Running efficient HIMA with 'MCP' penalty...")
      if (d$type != "continuous" || mediator.type != "gaussian") {
        stop("Efficient HIMA is only applicable to mediators and outcomes that are BOTH continuous and normally distributed.")
      }
      res <- hima_efficient(
        X = d$X,
        M = data.M,
        Y = d$Y,
        COV = d$COV,
        topN = NULL,
        scale = scale,
        FDRcut = sigcut,
        verbose = verbose,
        parallel = parallel,
        ncore = ncore
      )
      results <- .res_prep(res, method_text = "DACT method with BH-FDR", sigcut = sigcut)
    }
    
    if (quantile) {
      if (verbose) message("Running quantile HIMA with ", penalty, " penalty...")
      if (d$type != "continuous" || mediator.type != "gaussian") {
        stop("Quantile HIMA is only applicable to mediators and outcomes that are BOTH continuous and normally distributed.")
      }
      res <- hima_quantile(
        X = d$X,
        M = data.M,
        Y = d$Y,
        COV = d$COV,
        penalty = penalty,
        topN = NULL,
        scale = scale,
        Bonfcut = sigcut,
        verbose = verbose,
        parallel = parallel,
        ncore = ncore,
        ...
      )
      results <- .res_prep(res, method_text = "Bonferroni-adjusted p", sigcut = sigcut, q = TRUE)
    }
    
    if (longitudinal) {
      if (verbose) message("Running longitudinal survival HIMA with L1-penalized Cox lasso...")
      res <- hima_survival_long(
        X = d$X,
        M = data.M,
        tstart = d$Y$tstart,
        tstop = d$Y$tstop,
        status = d$Y$status,
        id = id_vec,
        COV = d$COV,
        topN = NULL,
        Bonfcut = sigcut,
        verbose = verbose,
        parallel = parallel,
        ncore = ncore
      )
      results <- .res_prep(res, method_text = "Bonferroni-adjusted p", sigcut = sigcut)
    }
  } else { # If not hima_efficient AND not hima_quantile AND not hima_survival_long
    
    if (penalty == "DBlasso") {
      if (d$type == "continuous") {
        if (mediator.type == "gaussian") {
          if (verbose) message("Running linear HIMA with de-biased Lasso penalty...")
          res <- hima_dblasso(
            X = d$X,
            M = data.M,
            Y = d$Y,
            COV = d$COV,
            topN = NULL,
            scale = scale,
            FDRcut = sigcut,
            verbose = verbose,
            parallel = parallel,
            ncore = ncore
          )
          results <- .res_prep(res, method_text = "HDMT pointwise FDR", sigcut = sigcut)
        } else if (mediator.type == "compositional") {
          if (verbose) message("Running compositional HIMA with de-biased Lasso penalty...")
          res <- hima_microbiome(
            X = d$X,
            OTU = data.M,
            Y = d$Y,
            COV = d$COV,
            FDRcut = sigcut,
            verbose = verbose,
            parallel = parallel,
            ncore = ncore
          )
          results <- .res_prep(res, method_text = "Hommel FDR", sigcut = sigcut)
        }
      } else if (d$type == "survival") {
        if (verbose) message("Running survival HIMA with de-biased Lasso penalty...")
        res <- hima_survival(
          X = d$X,
          M = data.M,
          OT = d$Y$OT,
          status = d$Y$status,
          COV = d$COV,
          topN = NULL,
          scale = scale,
          FDRcut = sigcut,
          verbose = verbose,
          parallel = parallel,
          ncore = ncore
        )
      }
      results <- .res_prep(res, method_text = "HDMT pointwise FDR", sigcut = sigcut)
    } else { # If penalty is not DBlasso
      res <- hima_classic(
        X = d$X,
        M = data.M,
        Y = d$Y,
        COV.XM = d$COV,
        Y.type = d$type,
        M.type = mediator.type,
        penalty = penalty,
        topN = NULL,
        scale = scale,
        Bonfcut = sigcut,
        verbose = verbose,
        parallel = parallel,
        ncore = ncore
      )
      results <- .res_prep(res, method_text = "Bonferroni-adjusted p", sigcut = sigcut)
    }
  }
  
  if (is.null(results) && verbose) message("No significant mediator found!")
  
  return(structure(results, class = "hima"))
}

#' @export
summary.hima <- function(object, desc = FALSE, ...) {
  if (is.null(object) || length(object$ID) == 0) {
    cat("No significant mediators identified.\n")
    return(invisible(NULL))
  }
  
  required_fields <- c("ID", "alpha", "beta", "alpha*beta", "rimp", "p-value")
  for (field in required_fields) {
    if (!field %in% names(object)) {
      stop(sprintf("The required field '%s' is missing from the object.", field))
    }
    if (!is.vector(object[[field]]) || is.null(object[[field]])) {
      stop(sprintf("The field '%s' must be a valid vector.", field))
    }
  }
  
  has_tau <- "tau" %in% names(object)
  
  if (has_tau) {
    result <- data.frame(
      ID = object$ID,
      alpha = object$alpha,
      beta = object$beta,
      `alpha*beta` = object$`alpha*beta`,
      `rimp` = object$rimp,
      `p-value` = object$`p-value`,
      tau = object$tau,  
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    result <- data.frame(
      ID = object$ID,
      alpha = object$alpha,
      beta = object$beta,
      `alpha*beta` = object$`alpha*beta`,
      `rimp` = object$rimp,
      `p-value` = object$`p-value`,
      check.names = FALSE,  
      stringsAsFactors = FALSE
    )
  }
  
  if (nrow(result) == 0) {
    cat("No significant mediators identified.\n")
    return(invisible(NULL))
  }
  
  cat("Summary of HIMA results:\n")
  cat("-----------------------------------\n")
  cat("Number of significant mediators:", nrow(result), "\n")
  
  cat("\nTop mediators (sorted by p-value):\n")
  sorted_result <- result[order(result$`p-value`), ]
  print(sorted_result)
  
  if (desc && !is.null(attr(object, "variable.labels"))) {
    cat("\nVariable Descriptions:\n")
    cat(paste(attr(object, "variable.labels"), collapse = "\n"))
  }
  
  return(invisible(NULL))
}



### Internal function: Result prepare function
.res_prep <- function(res, method_text, sigcut, q = FALSE) {
  if (is.null(res)) {
    return(NULL)
  }
  results <- data.frame(
    ID = res$Index,
    alpha = res$alpha_hat,
    beta = res$beta_hat,
    `alpha*beta` = res$IDE,
    `rimp` = res$rimp,
    `p-value` = res$pmax,
    check.names = FALSE
  )
  results <- results[order(results$`p-value`), ]
  attr_test <- c(
    "ID: Mediator ID",
    "alpha: Effect of exposure on mediator",
    "beta: Effect of mediator on outcome",
    "alpha*beta: Mediation (indirect) effect",
    "rimp: Relative importance of the mediator \n  out of all significant mediators",
    paste0("p-value: Joint raw p-value of significant mediator \n  selected based on ", method_text, " < ", sigcut)
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



### Internal function: recognize character variables and convert to dummy (only in hima)
.convert_to_dummies <- function(df) {
  char_cols <- sapply(df, is.character)
  df[char_cols] <- lapply(df[char_cols], as.factor)
  df_dummies <- as.data.frame(model.matrix(~ . - 1, data = df, na.action = na.pass))
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
  is_surv <- identical(as.character(f[[2]])[1], "Surv")
  
  if (is_surv) {
    response_vars <- all.vars(f[[2]])
    n_resp <- length(response_vars)
    if (!n_resp %in% c(2, 3)) {
      stop("Surv() must be specified as Surv(time, status) or Surv(tstart, tstop, status).")
    }
    
    if (n_resp == 2) {
      Y <- list(
        OT = data.pheno[[response_vars[1]]],
        status = data.pheno[[response_vars[2]]],
        tstart = NULL,
        tstop = NULL
      )
    } else {
      Y <- list(
        status = data.pheno[[response_vars[3]]],
        tstart = data.pheno[[response_vars[1]]],
        tstop = data.pheno[[response_vars[2]]]
      )
    }
    ind_vars <- all.vars(f[[3]])
    type <- "survival"
  } else {
    response_var <- as.character(f[[2]])
    ind_vars <- all.vars(f[[3]])
    Y <- data.pheno[[response_var]]
    type <- .data_type(Y)
  }
  if (!length(ind_vars)) {
    stop("No exposure/covariate specified on the right-hand side of the formula.")
  }
  
  X <- data.pheno[[ind_vars[1]]]
  
  if (length(ind_vars) > 1) {
    COV <- .convert_to_dummies(data.pheno[, ind_vars[-1], drop = FALSE])
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



### hima2 function has been renamed to hima with parameters updated.
hima2 <- function(formula,
                  data.pheno,
                  data.M,
                  outcome.family = c("gaussian", "binomial", "survival", "quantile"),
                  mediator.family = c("gaussian", "negbin", "compositional"),
                  penalty = c("DBlasso", "MCP", "SCAD", "lasso"),
                  efficient = FALSE,
                  scale = TRUE,
                  sigcut = 0.05,
                  verbose = FALSE,
                  ...) {
  .Deprecated(msg = "hima2() is deprecated. Please use hima() with updated parameters.")
  
  outcome.family <- match.arg(outcome.family)
  if (outcome.family == "quantile") quantile <- TRUE else quantile <- FALSE
  if (!missing(outcome.family)) {
    warning("'outcome.family' is deprecated. The outcome data type can now be automatically recognized.")
  }
  
  mediator.type <- match.arg(mediator.family)
  penalty <- match.arg(penalty)
  
  hima(formula,
       data.pheno,
       data.M,
       mediator.type = mediator.type,
       penalty = penalty,
       quantile = quantile,
       efficient = efficient,
       scale = scale,
       sigcut = sigcut,
       verbose = verbose,
       ...
  )
}
