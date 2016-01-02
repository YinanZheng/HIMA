#' Sure Independent Screening
#' 
#' \code{sis} is used to conduct sure independent screening using linear regression across high-dimensional variables as outcomes.
#' 
#' @param Y a data frame or matrix of genomic dataset. Rows represent samples, columns represent variables.
#' @param x.vect a vector of independent variable.
#' @param COV a data frame or matrix of covariates dataset. Default = \code{NULL}.
#' @param parallel logical. Enable parallel computing feature? Default = \code{TRUE}.
#' @param ncore number of cores to run parallel. Effective when paralle = \code{TRUE}. Default = 2.
#' 
#' @return a vector of estimates of mediators in their original order.
#' 
#' @seealso see \code{\link{hima}} using the results from \code{\link{sis}} as input to run PGS.
#'
#' @examples
#'
#' n <- 100  # sample size
#' p <- 5000 # the dimension of covariates
#' 
#' alpha <- rep(0, p) # the regression coefficients alpha (exposure --> mediators)
#' beta <- rep(0, p) # the regression coefficients beta (mediators --> outcome)
# 
#' alpha[1:3]  <- c(0.5, 0.5, 0.3) 
#' beta[1:3] <- c(0.5, 1.2, 0.3) # the first three makers are true mediators.
#' intercept <- runif(p,0,2)
#' 
#' set.seed(2016)
#' simdat = simHIMA(n, p, alpha, beta, intercept)

#' @export
sis <- function(Y, x.vect, COV=NULL, parallel=TRUE, ncore = 2)
{
  L.Y = ncol(Y); Y.names = colnames(Y)
  res.COV = NULL
  
  if (is.null(COV))
  { 
    datarun = data.frame(Yone = NA, x.vect=x.vect)
    modelstatement_MLR = modelstatement_GEE = Yone ~ x.vect
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~.,COV))[,-1]
    conf.names = colnames(COV)          
    datarun = data.frame(Yone = NA, x.vect=x.vect, COV = COV)
    modelstatement_MLR = eval(parse(text=(paste0("Yone ~ x.vect +", paste0(paste0("COV.",conf.names),collapse = "+")))))
  }
  
  modelstatement = modelstatement_MLR
  doOne <- eval(parse(text = doOneGen("try(lm(modelstatement, data = datarun))","1")))
    
  checkParallel("Sure Independent Screening", parallel, ncore, verbose = F)
  
  results <- foreach(n = idiv(L.Y, chunks = ncore), Y_chunk = iblkcol_lag(Y, chunks = ncore),.combine = 'c') %dopar% {
    do.call('rbind',lapply( seq_len(n), doOne, datarun, Y_chunk) )
  }
  
  names(results) = Y.names
  return(results)
}