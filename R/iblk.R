## Internal function: create iterator for bulk matrix by column

iblkcol_lag <- function(M, ...) {
  i <- 1
  it <- idiv(ncol(M), ...)
  
  nextEl <- function() {
    n <- nextElem(it)
    r <- seq(i, length=n)
    i <<- i + n
    M[,r, drop=FALSE]
  }
  obj <- list(nextElem=nextEl)
  class(obj) <- c('abstractiter', 'iter')
  obj
}