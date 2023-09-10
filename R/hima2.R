#' Advanced High-dimensional Mediation Analysis
#' 
#' \code{hima2} is an upgraded version of hima for estimating and testing high-dimensional mediation effects.
#' 
#' @param formula an object of class \code{formula}: a symbolic description of the overall effect model, i.e., 
#' \code{outcome ~ exposure + covariates}, to be fitted. Make sure the "exposure" is the variable of interest, which
#' must be listed as the first variable in the right hand side of the formula.
#' independent variable in the formula. The same covariates will be used in screening and penalized regression.
#' @param data.pheno a data frame containing exposure and covariates that are listed in the right hand side of the \code{formula}.
#' The variable names must match those listed in \code{formula}. By default \code{hima2} will scale \code{data.pheno}.
#' @param data.M a \code{data.frame} or \code{matrix} of high-dimensional mediators. Rows represent samples, columns 
#' represent variables. By default \code{hima2} will scale \code{data.M}.
#' @param outcome.family either \code{'gaussian'} (default, for normally distributed continuous outcome), \code{'binomial'} 
#' (for binay outcome), \code{'survival'} (for time-to-event outcome), or \code{'quantile'} (for quantile mediation analysis)
#' @param mediator.family either \code{'gaussian'} (default, for continuous mediators), \code{'negbin'} (i.e., negative binomial, 
#' for RNA-seq data as mediators), or \code{'compositional'} (for microbiome data as mediators), depending on the data type of 
#' high-dimensional mediators (\code{data.M}).
#' @param penalty the penalty to be applied to the model. Either \code{'DBlasso'} (De-biased LASSO, default), 
#' \code{'MCP'}, \code{'SCAD'}, or \code{'lasso'}. Please note, survival HIMA and microbiome HIMA can be only performed with \code{'DBlasso'}; 
#' Quantile HIMA cannot be performed with \code{'DBlasso'}.
#' @param topN an integer specifying the number of top markers from sure independent screening. 
#' Default = \code{NULL}. If \code{NULL}, \code{topN} will be \code{ceiling(2 * n/log(n))}, where \code{n} is the sample size. 
#' If the sample size is greater than topN (pre-specified or calculated), all mediators will be included in the test 
#' (i.e. a low-dimensional scenario).
#' @param scale logical. Should the function scale the data (exposure, mediators, and covariates)? Default = \code{TRUE}.
#' @param verbose logical. Should the function be verbose and shows the progression? Default = \code{FALSE}.
#' @param ... other arguments.
#' 
#' @return A data.frame containing mediation testing results of selected mediators. 
#' 
#' @references Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J, Just A, Colicino E, Vokonas P, Zhao L, 
#' Lv J, Baccarelli A, Hou L, Liu L. Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies. 
#' Bioinformatics. 2016. DOI: 10.1093/bioinformatics/btw351. PMID: 27357171. PMCID: PMC5048064
#' 
#' Perera C, Zhang H, Zheng Y, Hou L, Qu A, Zheng C, Xie K, Liu L. 
#' HIMA2: high-dimensional mediation analysis and its application in epigenome-wide DNA methylation data. 
#' BMC Bioinformatics. 2022. DOI: 10.1186/s12859-022-04748-1. PMID: 35879655. PMCID: PMC9310002
#' 
#' Zhang H, Zheng Y, Hou L, Zheng C, Liu L. Mediation Analysis for Survival Data with High-Dimensional Mediators. 
#' Bioinformatics. 2021. DOI: 10.1093/bioinformatics/btab564. PMID: 34343267. PMCID: PMC8570823
#' 
#' Zhang H, Chen J, Feng Y, Wang C, Li H, Liu L. Mediation effect selection in high-dimensional and compositional microbiome data. 
#' Stat Med. 2021. DOI: 10.1002/sim.8808. PMID: 33205470; PMCID: PMC7855955
#' 
#' Zhang H, Chen J, Li Z, Liu L. Testing for mediation effect with application to human microbiome data. 
#' Stat Biosci. 2021. DOI: 10.1007/s12561-019-09253-3. PMID: 34093887; PMCID: PMC8177450
#' 
#' Zhang H, Hong X, Zheng Y, Hou L, Zheng C, Wang X, Liu L. High-Dimensional Quantile Mediation Analysis with Application to a Birth 
#' Cohort Study of Mother–Newborn Pairs. Bioinformatics. 2023. (In press)
#' 
#' @examples
#' \dontrun{
#' # Note: In the following examples, M1, M2, and M3 are true mediators.
#' data(himaDat)
#' 
#' # Example 1 (continous outcome): 
#' head(himaDat$Example1$PhenoData)
#' 
#' e1 <- hima2(Outcome ~ Treatment + Sex + Age, 
#'       data.pheno = himaDat$Example1$PhenoData, 
#'       data.M = himaDat$Example1$Mediator,
#'       outcome.family = "gaussian",
#'       mediator.family = "gaussian",
#'       penalty = "DBlasso",
#'       scale = FALSE) # Disabled only for example data
#' e1
#' attributes(e1)$variable.labels
#' 
#' # Example 2 (binary outcome): 
#' head(himaDat$Example2$PhenoData)
#' 
#' e2 <- hima2(Disease ~ Treatment + Sex + Age, 
#'       data.pheno = himaDat$Example2$PhenoData, 
#'       data.M = himaDat$Example2$Mediator,
#'       outcome.family = "binomial",
#'       mediator.family = "gaussian",
#'       penalty = "DBlasso",
#'       scale = FALSE) # Disabled only for example data
#' e2
#' attributes(e2)$variable.labels
#' 
#' # Example 3 (time-to-event outcome): 
#' head(himaDat$Example3$PhenoData)
#' 
#' e3 <- hima2(Surv(Status, Time) ~ Treatment + Sex + Age, 
#'       data.pheno = himaDat$Example3$PhenoData, 
#'       data.M = himaDat$Example3$Mediator,
#'       outcome.family = "survival",
#'       mediator.family = "gaussian",
#'       penalty = "DBlasso",
#'       scale = FALSE) # Disabled only for example data
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
#'       penalty = "DBlasso",
#'       scale = FALSE) # Disabled only for example data
#' e4
#' attributes(e4)$variable.labels
#' 
#' #' # Example 5 (quantile mediation anlaysis): 
#' head(himaDat$Example5$PhenoData)
#' 
#' # Note that the function will prompt input for quantile level.
#' e5 <- hima2(Outcome ~ Treatment + Sex + Age, 
#'       data.pheno = himaDat$Example5$PhenoData, 
#'       data.M = himaDat$Example5$Mediator,
#'       outcome.family = "quantile",
#'       mediator.family = "gaussian",
#'       penalty = "MCP", # Quantile HIMA does not support DBlasso
#'       scale = FALSE, # Disabled only for example data
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
                  topN = NULL, 
                  scale = TRUE,
                  verbose = FALSE,
                  ...) 
{
  outcome.family <- match.arg(outcome.family)
  mediator.family <- match.arg(mediator.family)
  penalty <- match.arg(penalty)
  
  # Penalty check
  if (penalty == "DBlasso" & outcome.family == "quantile")
  {
    message("Note: Quantile HIMA does not support De-biased Lasso penalty. Switing to 'MCP' (default) ...")
    penalty = "MCP"
  }
  
  if (penalty != "DBlasso" & (outcome.family == "survival" | mediator.family == "compositional"))
  {
    message("Note: Survival HIMA and Compositional HIMA can be only performed using De-biased Lasso. Switing to 'DBlasso' ...")
    penalty = "DBlasso"
  }
  
  # DBlasso
  if (penalty == "DBlasso")
  {
    if (outcome.family %in% c("gaussian", "binomial"))
    {
      if(mediator.family %in% c("gaussian", "negbin"))
      {
        response_var <- as.character(formula[[2]]) 
        ind_vars <- all.vars(formula)[-1]
        
        Y <- data.pheno[,response_var]
        X <- data.pheno[,ind_vars[1]]
        
        if(length(ind_vars) > 1)
          COV <- data.pheno[,ind_vars[-1]] else COV <- NULL
        
        results <- dblassoHIMA(X = X, Y = Y, M = data.M, Z = COV, 
                               Y.family = outcome.family, 
                               topN = topN,
                               scale = scale, verbose = verbose)
        
        attr(results, "variable.labels") <- c("alpha: Effect of exposure on mediator", 
                                              "beta: Effect of mediator on outcome",
                                              "gamma: Total effect of exposure on outcome",
                                              "alpha*beta: Mediation effect",
                                              "% total effect: Percent of mediation effect out of the total effect",
                                              "p.joint: Joint raw p-value of selected significant mediator (based on FDR)")
      } else if (mediator.family == "compositional") {
        response_var <- as.character(formula[[2]]) 
        ind_vars <- all.vars(formula)[-1]
        
        Y <- data.pheno[,response_var]
        X <- data.pheno[,ind_vars[1]]
        
        if(length(ind_vars) > 1)
          COV <- data.pheno[,ind_vars[-1]] else COV <- NULL
        
        res <- microHIMA(X = X, Y = Y, OTU = data.M, COV = COV, FDRcut = 0.05, scale)
        results <- data.frame(alpha = res$alpha, alpha_se = res$alpha_se, 
                              beta = res$beta, beta_se = res$beta_se,
                              FDR = res$FDR, check.names = FALSE)
        rownames(results) <- res$ID
        attr(results, "variable.labels") <- c("alpha: Effect of exposure on mediator", 
                                              "alpha_se: Standard error of the effect of exposure on mediator",
                                              "beta: Effect of mediator on outcome",
                                              "beta_se: Standard error of the effect of mediator on outcome",
                                              "FDR: Hommel's false discovery rate")
      }
    } else if (outcome.family == "survival") {
      response_vars <- as.character(formula[[2]])[c(2,3)]
      ind_vars <- all.vars(formula)[-c(1,2)]
      
      X <- data.pheno[,ind_vars[1]]
      status <- data.pheno[, response_vars[1]]
      OT <- data.pheno[, response_vars[2]]
      
      if(length(ind_vars) > 1)
        COV <- data.pheno[,ind_vars[-1]] else COV <- NULL
      
      res <- survHIMA(X, COV, data.M, OT, status, FDRcut = 0.05, scale, verbose)
      
      results <- data.frame(alpha = res$alpha, alpha_se = res$alpha_se, 
                            beta = res$beta, beta_se = res$beta_se,
                            p.joint = res$p.joint, check.names = FALSE)
      rownames(results) <- res$ID
      attr(results, "variable.labels") <- c("alpha: Effect of exposure on mediator", 
                                            "alpha_se: Standard error of the effect of exposure on mediator",
                                            "beta: Effect of mediator on outcome",
                                            "beta_se: Standard error of the effect of mediator on outcome",
                                            "p.joint: Joint raw p-value of selected significant mediator (based on FDR)")
    }
  } else { # If penalty is not DBlasso
    if (outcome.family %in% c("gaussian", "binomial"))
    {
      response_var <- as.character(formula[[2]]) 
      ind_vars <- all.vars(formula)[-1]
      
      Y <- data.pheno[,response_var]
      X <- data.pheno[,ind_vars[1]]
      
      if(length(ind_vars) > 1)
        COV <- data.pheno[,ind_vars[-1]] else COV <- NULL
      
      results <- hima(X = X, Y = Y, M = data.M, COV.XM = COV, 
                      Y.family = outcome.family, M.family = mediator.family, 
                      penalty = penalty, topN = topN,
                      parallel = FALSE, ncore = 1, scale = scale, verbose = verbose)
      
      attr(results, "variable.labels") <- c("alpha: Effect of exposure on mediator", 
                                            "beta: Effect of mediator on outcome",
                                            "gamma: Total effect of exposure on outcome",
                                            "alpha*beta: Mediation effect",
                                            "% total effect: Percent of mediation effect out of the total effect",
                                            "Bonferroni.p: Bonferroni adjusted p value",
                                            "BH.FDR: Benjamini-Hochberg False Discovery Rate")
    } else if (outcome.family == "quantile") {
      # tau <- readline(prompt = "Enter quantile level(s) (between 0-1, multiple values accepted): ")
      # tau <- eval(parse(text = paste0("c(", tau, ")")))
      
      response_var <- as.character(formula[[2]]) 
      ind_vars <- all.vars(formula)[-1]
      
      Y <- data.pheno[,response_var]
      X <- data.pheno[,ind_vars[1]]
      
      if(length(ind_vars) > 1)
        COV <- data.pheno[,ind_vars[-1]] else COV <- NULL
      
      res <- qHIMA(X = X, M = data.M, Y = Y, Z = COV,
                   Bonfcut = 0.05, penalty = penalty, scale = scale, verbose = verbose, ...)
      
      results <- data.frame(alpha = res$alpha, alpha_se = res$alpha_se, 
                            beta = res$beta, beta_se = res$beta_se,
                            Bonferroni.p = res$Bonferroni.p, tau = res$tau, 
                            check.names = FALSE)
      rownames(results) <- paste0(res$ID, "-q", res$tau*100) 
      
      attr(results, "variable.labels") <- c("alpha: Effect of exposure on mediator", 
                                            "alpha_se: Standard error of the effect of exposure on mediator",
                                            "beta: Effect of mediator on outcome",
                                            "beta_se: Standard error of the effect of mediator on outcome",
                                            "Bonferroni.p: Bonferroni adjusted p value",
                                            "tau: Quantile level of the outcome")
    }
  }
  return(results)
}
