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
#' @param FDPcut FDP (false discovery proportions) cutoff applied to define and select significant mediators. Default = \code{0.05}. 
#' 
#' @return A data.frame containing mediation testing results of selected mediators (FDP < \code{FDPcut}). 
#' \itemize{
#'     \item{ID: }{index of selected significant mediator.}
#'     \item{alpha: }{coefficient estimates of exposure (X) --> mediators (M).}
#'     \item{alpha_se: }{standard error for alpha.}
#'     \item{beta: }{coefficient estimates of mediators (M) --> outcome (Y) (adjusted for exposure).}
#'     \item{beta_se: }{standard error for beta}
#'     \item{p_FDP: }{false discovery proportions of selected significant mediator.}
#' }
#' 
#' @references Zhang H, Chen J, Feng Y, Wang C, Li H, Liu L. Mediation effect selection in high-dimensional and compositional microbiome data. 
#' Stat Med. 2021. DOI: 10.1002/sim.8808. PMID: 33205470; PMCID: PMC7855955.
#' 
#' Zhang H, Chen J, Li Z, Liu L. Testing for mediation effect with application to human microbiome data. 
#' Stat Biosci. 2021. DOI: 10.1007/s12561-019-09253-3. PMID: 34093887; PMCID: PMC8177450.
#' 
#' @examples
#' \dontrun{
#' data(Example4)
#' head(Example4$PhenoData)
#' 
#' microHIMA.fit <- microHIMA(X = Example4$PhenoData$Treatment, 
#'                            Y = Example4$PhenoData$Outcome, 
#'                            OTU = Example4$Mediator, 
#'                            COV = Example4$PhenoData[, c("Sex", "Age")],
#'                            scale = FALSE)
#' microHIMA.fit
#' }
#' 
#' @export
microHIMA <- function(X, Y, OTU, COV = NULL, FDPcut = 0.05, scale = TRUE){
  
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
  
  ## The FDP method
  set <- which(P_adj_DLASSO < FDPcut)
  hom <- hommel::hommel(P_adj_DLASSO, simes = FALSE)
  N1 <- hommel::discoveries(hom, set, incremental = TRUE, alpha=0.05)
  
  if (length(set) > 0){
    L <- length(set)
    N2 <- matrix(0,1,L)
    N2[2:L] <- N1[1:(L-1)]
  }
  
  N0 <- N1 - N2
  
  ID_FDR <- set[which(N0 > 0)]

  out_result <- data.frame(ID = M_ID_name[ID_FDR], 
                           alpha = alpha_EST[ID_FDR], 
                           alpha_se = alpha_SE[ID_FDR], 
                           beta = beta_EST[ID_FDR], 
                           beta_se = beta_SE[ID_FDR],
                           p_FDP = P_adj_DLASSO[ID_FDR])
  
  message("Done!", "     (", format(Sys.time(), "%X"), ")")
  
  return(out_result)
}
