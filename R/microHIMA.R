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
#' ## Generate simulated survival data
#' n <- 200 # Sample size
#' p <- 25 # Number of microbiome
#' Treatment = rbinom(n, 1, 0.2) # binary outcome
#'
#' ## Generate two covariates, one binary, one continuous
#' covariates = cbind(sample(c(1,0), n, replace = TRUE), rnorm(n))
#' 
#' ## parameters
#' beta0 = as.numeric(matrix(0, 1, p))
#' betaT = rep(0, p)
#' betaT[c(1, 2, 3)] = c(1, 1.2, 1.5) # let the first three are non-zero
#' betaX = matrix(0, p, 2)
#' 
#' alpha0 = 0
#' alphaT = 1
#' alphaZ = alphaC = rep(0, p)
#' alphaZ[c(1, 2, 3)] = c(1.3, -0.7, -0.6) # let the first three are non-zero for response
#' alphaX = c(0.5, 0.5)
#' 
#' ## Generate microbiome data
#' X = cbind(rep(1, n), covariates, Treatment) # n * (1 + q + p)
#' b = cbind(beta0, betaX, betaT) # p * (1 + q + p)
#' gamma.simu = exp(X %*% t(b)) # n * p
#' otu.com = t(apply(gamma.simu, 1, HIMA:::rdirichlet, n = 1)) # Dirichlet distribution
#' 
#' ## Generate outcome data
#' X = cbind(rep(1, n),Treatment, covariates, log(otu.com), log(otu.com) * Treatment)
#' b = c(alpha0, alphaT, alphaX, alphaZ, alphaC)
#' outcome = c(b %*% t(X) + rnorm(n, mean = 0, sd = 1))
#' exposure = t(t(Treatment))
#' 
#' \dontrun{
#' microHIMA.fit <- microHIMA(X = exposure, Y = outcome, OTU = otu.com, COV = covariates)
#' microHIMA.fit
#' }
#' 
#' @export
microHIMA <- function(X, Y, OTU, COV = NULL, FDPcut = 0.05){
  
  M_raw <- OTU
  X <- cbind(X, COV)
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
  
  message("Running debiased Lasso...", "     (", Sys.time(), ")")
  
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
    
    MT <- matrix(as.numeric(scale(MT)), nrow(scale(MT)), ncol(scale(MT)))
    X <- matrix(as.numeric(scale(X)), nrow(scale(X)), ncol(scale(X)))
    MX <- cbind(MT, X)

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

  out_result <- data.frame(ID = ID_FDR, 
                           alpha = alpha_EST[ID_FDR], 
                           alpha_se = alpha_SE[ID_FDR], 
                           beta = beta_EST[ID_FDR], 
                           beta_se = beta_SE[ID_FDR],
                           p_FDP = P_adj_DLASSO[ID_FDR])
  
  message("Done!", "     (", Sys.time(), ")")
  
  return(out_result)
}
