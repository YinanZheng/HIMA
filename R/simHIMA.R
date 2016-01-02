#' Simulation Data Generator for High-dimensional Mediation Analyais
#' 
#' \code{simHIMA} is used to simulate data for High-dimensional Mediation Analyais
#' 
#' @param n an integer specifying sample size.
#' @param p an integer specifying the dimension of mediators.
#' @param alpha a numeric vector specifying the regression coefficients alpha (exposure --> mediators).
#' @param beta a numeric vector specifying the regression coefficients beta (mediators --> outcome).
#' @param seed an integer specifying a seed for random number generation.
#' 
#' @seealso see \code{\link{hima}} to run HMA algorithm.
#' 
#' @examples
#' n <- 100  # sample size
#' p <- 5000 # the dimension of covariates
#' 
#' alpha <- rep(0, p) # the regression coefficients alpha (exposure --> mediators)
#' beta <- rep(0, p) # the regression coefficients beta (mediators --> outcome)
# 
#' alpha[1:3] <- c(0.5, 0.5, 0.3) 
#' beta[1:3] <- c(0.5, 1.2, 0.3) # the first three makers are true mediators.
#' 
#' simdat = simHMA(n, p, alpha, beta, seed=2016) # Generate simulation data
#' 
#' @export
simHIMA <- function(n, p, alpha, beta, seed)
{
  set.seed(seed)
  intercept <- runif(p,0,2)
  M <- matrix(0,n,p) #  high-dimensional mediator matrix with n by p
  colnames(M) = paste0("M",1:p)
  x <- rnorm(n, mean=1, sd=1.5)
  e <- matrix(0,1,n)
  for (i in 1:p){
    e <- rnorm(n, 0, 1.2)
    y <- alpha[i]*x + intercept[i] + e
    M[,i] <- t(t(y))
  }
  X <- t(t(x)) #  independent confounder X with n by 1
  X_M <- cbind(M,X) # matrix combine M with X, 
  E <- rnorm(n,0,1)
  B <- c(beta, 0.5) # combind coefficient for X = 0.5 to beta
  Y <-  1 + X_M%*%t(t(B)) + t(t(E)) #  the response with  n by 1
  return(list(Y = Y, M = M, X = X, n = n, p = p))
}