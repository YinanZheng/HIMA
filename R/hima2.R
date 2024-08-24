#' High-dimensional Mediation Analysis
#' 
#' \code{hima2} is a wrapper function designed to perform various HIMA methods for estimating and testing high-dimensional mediation effects.
#' \code{hima2} can automatically select the appropriate HIMA method based on the outcome and mediator data type specified by the user.
#' 
#' @param formula an object of class \code{formula}: a symbolic description of the overall effect model, i.e., 
#' \code{outcome ~ exposure + covariates}, to be fitted. Make sure the "exposure" is the variable of interest, which
#' must be listed as the first variable in the right hand side of the formula.
#' @param data.pheno a data frame containing exposure and covariates that are listed in the right hand side of the \code{formula}.
#' The variable names must match those listed in \code{formula}. By default \code{hima2} will scale \code{data.pheno}.
#' @param data.M a \code{data.frame} or \code{matrix} of high-dimensional mediators. Rows represent samples, columns 
#' represent variables. By default \code{hima2} will scale \code{data.M}.
#' @param outcome.family either \code{'gaussian'} (default, for normally distributed continuous outcome), \code{'binomial'} 
#' (for binary outcome), \code{'survival'} (for time-to-event outcome), or \code{'quantile'} (for quantile mediation analysis)
#' @param mediator.family either \code{'gaussian'} (default, for continuous mediators), \code{'negbin'} (i.e., negative binomial, 
#' for RNA-seq data as mediators), or \code{'compositional'} (for microbiome data as mediators), depending on the data type of 
#' high-dimensional mediators (\code{data.M}).
#' @param penalty the penalty to be applied to the model. Either \code{'DBlasso'} (De-biased LASSO, default), 
#' \code{'MCP'}, \code{'SCAD'}, or \code{'lasso'}. Please note, survival HIMA and microbiome HIMA can be only performed with \code{'DBlasso'}; 
#' Quantile HIMA cannot be performed with \code{'DBlasso'}. Not applicable for efficient HIMA (when \code{efficient = TRUE}), as it will always apply \code{'MCP'}. 
#' @param efficient use efficient HIMA (\code{eHIMA}). Only applicable for linear HIMA with continuous outcome (i.e., \code{outcome.family = "gaussian"} 
#' and \code{mediator.family = "gaussian"}). Default = \code{FALSE}.
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
#' # Example 1 (continous outcome - linear hima): 
#' head(himaDat$Example1$PhenoData)
#' 
#' e1 <- hima2(Outcome ~ Treatment + Sex + Age, 
#'       data.pheno = himaDat$Example1$PhenoData, 
#'       data.M = himaDat$Example1$Mediator,
#'       outcome.family = "gaussian",
#'       mediator.family = "gaussian",
#'       penalty = "MCP", # Can be "DBlasso" for dblassoHIMA
#'       scale = FALSE) # Disabled only for simulation data
#' e1
#' attributes(e1)$variable.labels
#' 
#' # Efficient HIMA (only applicable to mediators and outcomes that are 
#' # both continuous and normally distributed.)
#' e1e <- hima2(Outcome ~ Treatment + Sex + Age,
#'              data.pheno = himaDat$Example1$PhenoData,
#'              data.M = himaDat$Example1$Mediator,
#'              outcome.family = "gaussian",
#'              mediator.family = "gaussian",
#'              efficient = TRUE,
#'              scale = FALSE) # Disabled only for simulation data
#' e1e
#' attributes(e1e)$variable.labels
#' 
#' # Example 2 (binary outcome - logistic hima): 
#' head(himaDat$Example2$PhenoData)
#' 
#' e2 <- hima2(Disease ~ Treatment + Sex + Age, 
#'       data.pheno = himaDat$Example2$PhenoData, 
#'       data.M = himaDat$Example2$Mediator,
#'       outcome.family = "binomial",
#'       mediator.family = "gaussian",
#'       penalty = "MCP",
#'       scale = FALSE) # Disabled only for simulation data
#' e2
#' attributes(e2)$variable.labels
#' 
#' # Example 3 (time-to-event outcome - survival hima): 
#' head(himaDat$Example3$PhenoData)
#' 
#' e3 <- hima2(Surv(Status, Time) ~ Treatment + Sex + Age, 
#'       data.pheno = himaDat$Example3$PhenoData, 
#'       data.M = himaDat$Example3$Mediator,
#'       outcome.family = "survival",
#'       mediator.family = "gaussian",
#'       penalty = "DBlasso",
#'       scale = FALSE) # Disabled only for simulation data
#' e3
#' attributes(e3)$variable.labels
#' 
#' # Example 4 (compositional data as mediator, e.g., microbiome): 
#' head(himaDat$Example4$PhenoData)
#' 
#' e4 <- hima2(Outcome ~ Treatment + Sex + Age, 
#'       data.pheno = himaDat$Example4$PhenoData, 
#'       data.M = himaDat$Example4$Mediator,
#'       outcome.family = "gaussian",
#'       mediator.family = "compositional",
#'       penalty = "DBlasso") # Scaling is always enabled for microHIMA
#' e4
#' attributes(e4)$variable.labels
#' 
#' #' # Example 5 (quantile mediation anlaysis - quantile hima): 
#' head(himaDat$Example5$PhenoData)
#' 
#' # Note that the function will prompt input for quantile level.
#' e5 <- hima2(Outcome ~ Treatment + Sex + Age, 
#'       data.pheno = himaDat$Example5$PhenoData, 
#'       data.M = himaDat$Example5$Mediator,
#'       outcome.family = "quantile",
#'       mediator.family = "gaussian",
#'       penalty = "MCP", # Quantile HIMA does not support DBlasso
#'       scale = FALSE, # Disabled only for simulation data
#'       tau = c(0.3, 0.5, 0.7)) # Specify multiple quantile level
#' e5
#' attributes(e5)$variable.labels
#' }
#'                   
#' @export
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
                  ...) 
{
  outcome.family <- match.arg(outcome.family)
  mediator.family <- match.arg(mediator.family)
  penalty <- match.arg(penalty)
  
  if (Sigcut > 1) Sigcut = 1
  
  if(is.null(colnames(data.M))) colnames(data.M) <- seq_len(ncol(data.M))
  
  # Penalty check
  if (penalty == "DBlasso")
  {
    if(outcome.family == "quantile") 
    {
      message("Note: Quantile HIMA does not support De-biased Lasso penalty. Switing to 'MCP' ...")
      penalty = "MCP"
    }
    if(outcome.family == "binomial") 
    {
      message("Note: Logistic HIMA does not support De-biased Lasso penalty. Switing to 'MCP' ...")
      penalty = "MCP"
    }
  }
    
  if (penalty != "DBlasso")
  {
    if(outcome.family == "survival")
    {
      message("Note: Survival HIMA can be only performed using De-biased Lasso. Switing to 'DBlasso' ...")
      penalty = "DBlasso"
    }
    if(mediator.family == "compositional")
    {
      message("Note: Compositional HIMA can be only performed using De-biased Lasso. Switing to 'DBlasso' ...")
      penalty = "DBlasso"
    }
  }
  
  ######
  
  # eHIMA
  if (efficient == TRUE)
  {
    if(outcome.family != "gaussian" | mediator.family != "gaussian") 
      stop("Efficient HIMA is only applicable to mediators and outcomes that are BOTH continuous and normally distributed.")
    
    d <- .f_extract(data.pheno = data.pheno, f = formula)
    
    res <- eHIMA(X = d$X, 
                 Y = d$Y, 
                 M = data.M, 
                 COV = d$COV, 
                 topN = NULL,
                 scale = scale,
                 FDRcut = Sigcut,
                 verbose = verbose)
    
    results <- .res_prep(res, method_text = "DACT method with BH-FDR", Sigcut = Sigcut)

  } else {
    # DBlasso
    if (penalty == "DBlasso")
    {
      if (outcome.family %in% c("gaussian", "binomial"))
      {
        if(mediator.family %in% c("gaussian", "negbin"))
        {
          d <- .f_extract(data.pheno = data.pheno, f = formula)
          
          res <- dblassoHIMA(X = d$X, 
                             Y = d$Y, 
                             M = data.M, 
                             COV = d$COV, 
                             Y.family = outcome.family, 
                             topN = NULL,
                             FDRcut = Sigcut,
                             scale = scale, 
                             verbose = verbose)
          
          results <- .res_prep(res, method_text = "HDMT pointwise FDR", Sigcut = Sigcut)
          
        } else if (mediator.family == "compositional") 
        {
          d <- .f_extract(data.pheno = data.pheno, f = formula)
          
          res <- microHIMA(X = d$X, 
                           Y = d$Y, 
                           OTU = data.M, 
                           COV = d$COV, 
                           FDRcut = Sigcut, 
                           verbose = verbose)
          
          results <- .res_prep(res, method_text = "Hommel FDR", Sigcut = Sigcut)
        }
      } else if (outcome.family == "survival") 
      {
        d <- .f_extract(data.pheno = data.pheno, f = formula, tte = TRUE)
        
        res <- survHIMA(X = d$X, 
                        M = data.M, 
                        COV = d$COV, 
                        OT = d$Y$OT, 
                        status = d$Y$status, 
                        FDRcut = Sigcut, 
                        scale = scale, 
                        verbose = verbose)
        
        results <- .res_prep(res, method_text = "HDMT pointwise FDR", Sigcut = Sigcut)
      }
    } else { # If penalty is not DBlasso
      if (outcome.family %in% c("gaussian", "binomial"))
      {
        d <- .f_extract(data.pheno = data.pheno, f = formula)
        
        res <- hima(X = d$X, 
                    Y = d$Y, 
                    M = data.M, 
                    COV.XM = d$COV, 
                    Y.family = outcome.family, 
                    M.family = mediator.family, 
                    penalty = penalty, 
                    topN = NULL,
                    parallel = FALSE, 
                    ncore = 1, 
                    scale = scale, 
                    Bonfcut = Sigcut,
                    verbose = verbose)
        
        results <- .res_prep(res, method_text = "Bonferroni-adjusted p", Sigcut = Sigcut)

      } else if (outcome.family == "quantile") 
      {
        # tau <- readline(prompt = "Enter quantile level(s) (between 0-1, multiple values accepted): ")
        # tau <- eval(parse(text = paste0("c(", tau, ")")))
        
        d <- .f_extract(data.pheno = data.pheno, f = formula)
        
        res <- qHIMA(X = d$X, 
                     M = data.M, 
                     Y = d$Y, 
                     COV = d$COV,
                     penalty = penalty, 
                     scale = scale, 
                     Bonfcut = 0.05,
                     verbose = verbose,
                     ...)
        
        results <- .res_prep(res, method_text = "Bonferroni-adjusted p", Sigcut = Sigcut, q = TRUE)
      }
    }
  }
  
  if(is.null(results)) message("No significant mediator found!")
  
  return(results)
}


### Internal function: Result prepare function
.res_prep <- function(res, method_text, Sigcut, q = FALSE)
{
  if(is.null(res)) return(NULL)
  results <- data.frame(ID = res$Index,
                        alpha = res$alpha_hat, 
                        beta = res$beta_hat, 
                        `alpha*beta` = res$IDE,
                        `Relative Importance (%)` = res$rimp,
                        `p-value` = res$pmax,
                        check.names = FALSE)
  results <- results[order(results$`p-value`), ]
  attr_test <- c("ID: Mediator ID",
                 "alpha: Effect of exposure on mediator", 
                 "beta: Effect of mediator on outcome",
                 "alpha*beta: Mediation (indirect) effect",
                 "Relative Importance (%): Relative importance of the mediator out of all significant mediators",
                 paste0("p-value: Joint raw p-value of significant mediator selected based on ", method_text, " < ", Sigcut))
  if(q) 
  {
    results <- data.frame(results, 
                          tau = res$tau, 
                          check.names = FALSE)
    attr_test <- c(attr_test,
                   "tau: Quantile level of the outcome")
  }
  attr(results, "variable.labels") <- attr_test
  rownames(results) <- NULL
  return(results)
}



### Internal function: recognize character variables and convert to dummy (only in hima2)
.convert_to_dummies <- function(df) {
  char_cols <- sapply(df, is.character)
  df[char_cols] <- lapply(df[char_cols], as.factor)
  df_dummies <- as.data.frame(model.matrix(~ . - 1, data = df))
  non_char_cols <- names(df)[!char_cols]
  df_dummies[non_char_cols] <- df[non_char_cols]
  return(df_dummies)
}



### Internal function: Formula extractor function
.f_extract <- function(data.pheno, f, tte = FALSE)
{
  if(tte)
  {
    response_vars <- as.character(f[[2]])[c(2,3)]
    ind_vars <- all.vars(f)[-c(1,2)]
    Y <- list(status = data.pheno[, response_vars[1]], 
              OT = data.pheno[, response_vars[2]])
  } else {
    response_var <- as.character(f[[2]]) 
    ind_vars <- all.vars(f)[-1]
    Y <- data.pheno[,response_var]
  }
  
  X <- data.pheno[,ind_vars[1]]
  
  if(length(ind_vars) > 1)
    COV <- .convert_to_dummies(data.pheno[,ind_vars[-1]]) else COV <- NULL
  
  return(list(Y = Y, 
              X = X, 
              COV = COV))
}
