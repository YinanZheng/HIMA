#' Simulation Data Generator for High-dimensional Mediation Analyais
#' 
#' \code{simHIMA} is used to generate simulation data for high-dimensional mediation analyais.
#' 
#' @param n an integer specifying sample size.
#' @param p an integer specifying the dimension of mediators.
#' @param alpha a numeric vector specifying the regression coefficients alpha (exposure --> mediators).
#' @param beta a numeric vector specifying the regression coefficients beta (mediators --> outcome).
#' @param binaryOutcome logical. Should the simulated outcome variable be binary?
#' @param seed an integer specifying a seed for random number generation.
#' 
#' @seealso see \code{\link{hima}} to run HIMA.
#' 
#' @examples
#' n <- 300  # sample size
#' p <- 10000 # the dimension of covariates
#' 
#' # the regression coefficients alpha (exposure --> mediators)
#' alpha <- rep(0, p) 
#' 
#' # the regression coefficients beta (mediators --> outcome)
#' beta <- rep(0, p) 
#' 
#' # the first four markers are true mediators.
#' alpha[1:4] <- c(0.45,0.5,0.55,0.6)
#' beta[1:4] <- c(0.5,0.45,0.4,0.35)
#' 
#' alpha[7:8] <- 0.5
#' beta[5:6] <- 0.5
#' 
#' # Generate simulation data
#' simdat = simHIMA(n, p, alpha, beta, seed=2016) 
#' 
#' @export
simHIMA <- function(n, p, alpha, beta, binaryOutcome = FALSE, seed) {
  set.seed(seed)
  ck <- t(runif(p, 0, 2))
  
  M <- matrix(0, n, p)
  X <- t(t(rnorm(n, mean = 0, sd = 1.5)))
  
  for (i in 1:n) {
    e <- rnorm(p, 0, 1.2)
    M[i, ] <- ck + X[i] * alpha + e  #  The mediators matrix n times p
  }
  
  colnames(M) <- paste0("M", 1:ncol(M))
  
  XM <- cbind(X, M)  #  [X M]
  B <- c(0.5, beta)  # (p+1) times 1
  E <- rnorm(n, 0, 1)
  Y <- 0.5 + XM %*% t(t(B)) + t(t(E))  #  the response  n times 1
  
  if(binaryOutcome)
    Y <- matrix(rbinom(n, 1, 1/(1+exp(-Y))), nrow = n)
  
  return(list(Y = Y, M = M, X = X, n = n, p = p))
}