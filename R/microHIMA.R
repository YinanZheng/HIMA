# This is the main function for our proposed method for high-dimensional compositional microbiome mediation analysis
#' High-dimensional mediation analysis for compositional microbiome data
#' 
#' \code{microHIMA} is used to estimate and test high-dimensional mediation effects for compositional microbiome data.
#' 
#' @param X a vector of exposure. 
#' @param Y a vector of outcome.
#' @param OTU a \code{data.frame} or \code{matrix} of high-dimensional compositional OTUs (mediators). Rows represent samples, 
#' columns represent variables.
#' @param COV a \code{data.frame} or \code{matrix} of adjusting covariates. Rows represent samples, columns represent microbiome variables. 
#' Can be \code{NULL}.
#' @param scale logical. Should the function scale the data? Default = \code{TRUE}.
#' @param FDRcut Hommel FDR cutoff applied to select significant mediators. Default = \code{0.05}. 
#' @param verbose logical. Should the function be verbose? Default = \code{FALSE}.
#' 
#' @return A data.frame containing mediation testing results of significant mediators (FDR <\code{FDRcut}). 
#' \describe{
#'     \item{Index: }{mediation name of selected significant mediator.}
#'     \item{alpha_hat: }{coefficient estimates of exposure (X) --> mediators (M).}
#'     \item{alpha_se: }{standard error for alpha.}
#'     \item{beta_hat: }{coefficient estimates of mediators (M) --> outcome (Y) (adjusted for exposure).}
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
#' data(himaDat)
#' 
#' head(himaDat$Example4$PhenoData)
#' 
#' microHIMA.fit <- microHIMA(X = himaDat$Example4$PhenoData$Treatment, 
#'                            Y = himaDat$Example4$PhenoData$Outcome, 
#'                            OTU = himaDat$Example4$Mediator, 
#'                            COV = himaDat$Example4$PhenoData[, c("Sex", "Age")],
#'                            scale = FALSE)
#' microHIMA.fit
#' }
#' 
#' @export
microHIMA <- function(X, Y, OTU, COV = NULL, 
                      FDRcut = 0.05, 
                      scale = TRUE, 
                      verbose = FALSE)
{
  
  X <- matrix(X, ncol = 1)
  
  M_raw <- as.matrix(OTU)
  
  M_ID_name <- colnames(M_raw)
  if(is.null(M_ID_name)) M_ID_name <- seq_len(ncol(M_raw))
  
  if(!is.null(COV))
  {COV <- as.matrix(COV); X <- cbind(X, COV)}
  
  Y <- Y - mean(Y)
  
  M <- M_raw
  n <- dim(M)[1]
  d <- dim(M)[2]
  Index_S <- matrix(0,1,d)
  P_b_raw <-  matrix(0,1,d)
  P_a_raw <-  matrix(0,1,d)
  alpha_EST <- matrix(0,1,d)
  alpha_SE  <- matrix(0,1,d)
  beta_EST <- matrix(0,1,d)
  beta_SE <-  matrix(0,1,d)
  P_raw_DLASSO <- matrix(0,1,d)
  M1 <- t(t(M_raw[,1]))
  
  message("Step 1: Isometric Log-ratio Transformation and De-biased Lasso estimates ...", "  (", format(Sys.time(), "%X"), ")")
  
  if(verbose)
  {
    if(is.null(COV)) 
    {message("        No covariate was adjusted.")} 
    else
    {message("        Adjusting for covariate(s): ", paste0(colnames(COV), collapse = ", "))}
  }
  
  for (k in 1:d){
    M <- M_raw
    M[,1] <- M[,k]
    M[,k] <- M1 
    MT <- matrix(0,n,d-1)
    for (i in 1:n){
      for (j in 1:(d-1)){
        C_1 <- sqrt((d-j)/(d-j+1))
        C_2 <- prod(M[i,(j+1):d]^(1/(d-j)))
        MT[i,j] <- C_1*log(M[i,j]/C_2)
      }
    }
    
    MT <- matrix(as.numeric(MT), nrow(MT))
    MX <- cbind(MT, X)
    
    if(scale) MX <- scale(MX)
    
    fit.dlasso  <- DLASSO_fun(MX, Y)
    
    beta_est <- fit.dlasso[1]
    beta_se  <- fit.dlasso[2]
    P_b <-  2*(1-pnorm(abs(beta_est/beta_se),0,1))
    beta_EST[k] <-  beta_est
    beta_SE[k]  <- beta_se
    
    lm.fit <- stats::lm(MT[,1]~X)
    lm.out <- summary(lm.fit)
    alpha_est <- lm.out$coefficients[2,1]
    alpha_se <- lm.out$coefficients[2,2]
    P_a <-  2*(1-pnorm(abs(alpha_est/alpha_se),0,1))
    P_raw_DLASSO[k] <- max(P_a,P_b)
    alpha_EST[k] <- alpha_est
    alpha_SE[k] <- alpha_se
    
  } #the end of k
  
  P_adj_DLASSO <- as.numeric(P_raw_DLASSO)
  
  message("Step 2: Closted testing-based procedure ...", "     (", format(Sys.time(), "%X"), ")")
  
  ## The FDR method
  set <- which(P_adj_DLASSO < FDRcut)
  hom <- hommel::hommel(P_adj_DLASSO, simes = FALSE)
  N1 <- hommel::discoveries(hom, set, incremental = TRUE, alpha=0.05)
  
  if (length(set) > 0){
    L <- length(set)
    N2 <- matrix(0,1,L)
    N2[2:L] <- N1[1:(L-1)]
  }
  
  N0 <- N1 - N2
  
  ID_FDR <- set[which(N0 > 0)]
  
  IDE <- alpha_EST[ID_FDR] * beta_EST[ID_FDR]
  
  if (length(ID_FDR) > 0){
    out_result <- data.frame(Index = M_ID_name[ID_FDR], 
                             alpha_hat = alpha_EST[ID_FDR], 
                             alpha_se = alpha_SE[ID_FDR], 
                             beta_hat = beta_EST[ID_FDR], 
                             beta_se = beta_SE[ID_FDR],
                             IDE = IDE, 
                             rimp = abs(IDE)/sum(abs(IDE)) * 100, 
                             pmax = P_adj_DLASSO[ID_FDR])
    if(verbose) message(paste0("        ", length(ID_FDR), " significant mediator(s) identified."))
  } else {
    if(verbose) message("        No significant mediator identified.")
    out_result = NULL
  }
  
  message("Done!", "     (", format(Sys.time(), "%X"), ")")
  
  return(out_result)
}
