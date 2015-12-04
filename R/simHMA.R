#' Simulation Data Generator for High-dimensional Mediation Analyais
#' 
#' \code{simHMA} is used to simulate data for High-dimensional Mediation Analyais
#' 
#' @param n an integer specifying sample size.
#' @param p an integer specifying the dimension of covariates.
#' @param alpha a numeric vector specifying the regression coefficients alpha
#' @param beta a numeric vector specifying the regression coefficients beta
#' 
#' @seealso see \code{\link{hma}} to run HMA algorithm.
#' 
#' @examples
#' n <- 100  # sample size
#' p <- 5000 # the dimension of covariates
#' 
#' alpha <- rep(0, p) # the regression coefficients alpha
#' beta <- rep(0, p) # the regression coefficients beta
#' 
#' alpha[1:3]  <- c(0.5, 0.5, 0.3)
#' beta[1:3] <- c(0.5, 1.2, 0.3)
#'
#' simdat = simHMA(n, p, alpha, beta)

simHMA <- function(n, p, a, beta)
{
  M <- matrix(0,n,p) #  high-dimensional mediator matrix with n by p
  x <- rnorm(n, mean=1, sd=1.5)
  e <- matrix(0,1,n)
  for (i in 1:p){
    e <- rnorm(n, 0, 1.2)
    y <- a[i]*x + c[i] + e
    M[,i] <- t(t(y))
  }
  X <- t(t(x)) #  independent confounder X with n by 1
  X_M <- cbind(M,X) # matrix combine M with X, 
  E <- rnorm(n,0,1)
  B <- c(b, 0.5) # combind coefficient for X = 0.5 to beta
  Y <-  1 + X_M%*%t(t(B)) + t(t(E)) #  the response with  n by 1
  return(list(Y = Y, M = M, X = X, n = n, p = p))
}