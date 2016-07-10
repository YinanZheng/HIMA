#' Simulation Data Generator for High-dimensional Mediation Analyais
#' 
#' \code{simHIMA} is used to generate simulation data for high-dimensional mediation analyais.
#' 
#' @param n an integer specifying sample size.
#' @param p an integer specifying the dimension of mediators.
#' @param alpha a numeric vector specifying the regression coefficients alpha (exposure --> mediators).
#' @param beta a numeric vector specifying the regression coefficients beta (mediators --> outcome).
#' @param seed an integer specifying a seed for random number generation.
#' 
#' @seealso see \code{\link{hima}} to run HIMA.
#' 
#' @examples
#' n <- 300  # sample size
#' p <- 10000 # the dimension of covariates
#' 
#' alpha <- rep(0, p) # the regression coefficients alpha (exposure --> mediators)
#' beta <- rep(0, p) # the regression coefficients beta (mediators --> outcome)
# 
#' alpha[1:4] <- c(0.45,0.55,0.5,0.55)
#' beta[1:4] <- c(0.40,0.55,0.5,0.40) # the first four makers are true mediators.
#' 
#' alpha[7:8] <- 0.55
#' beta[5:6] <- 0.50
#' 
#' simdat = simHIMA(n, p, alpha, beta, seed=2016) # Generate simulation data
#' 
#' @export
simHIMA <- function(n, p, alpha, beta, seed)
{
  set.seed(seed)
  ck <- t(runif(p,0,2))
  
  M <- matrix(0,n,p) 
  X  <-  t(t(rnorm(n, mean=0, sd=1.5)))
  
  for (i in 1:n){
    e <- rnorm(p, 0, 1.2)
    M[i,] <- ck + X[i]*alpha + e  #  The mediators matrix n times p
  }
 
  colnames(M) = paste0("M",1:ncol(M))
  
  XM <- cbind(X,M) #  [X M]
  B <- c(0.5,beta) # (p+1) times 1
  E <- rnorm(n,0,1)
  Y <-  0.5 + XM%*%t(t(B)) + t(t(E)) #  the response  n times 1

  return(list(Y = Y, M = M, X = X, n = n, p = p))
}