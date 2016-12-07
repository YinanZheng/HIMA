# Internal function: parallel computing check

checkParallel <- function(program.name, parallel, ncore, verbose = TRUE) {
  if (parallel & ncore > 1) {
    if (ncore > parallel::detectCores()) {
      message("You requested ", ncore, " cores. There are only ", 
              parallel::detectCores(), " in your machine!")
      ncore <- parallel::detectCores()
    }
    if (verbose) 
      message("    Running ", program.name, " with ", ncore, " cores in parallel...   (", 
              Sys.time(), ")")
      doParallel::registerDoParallel(ncore)
  } else {
    if (verbose) 
      message("    Running ", program.name, " with single core...   (", 
              Sys.time(), ")")
    registerDoSEQ()
  }
}

## Internal function: doOne code generater

doOneGen <- function(model.text, colind.text) {
  L <- length(eval(parse(text = colind.text)))
  script <- paste0("doOne <- function(i, datarun, Ydat){datarun$Mone <- Ydat[,i]; model <- ", 
                   model.text, ";if('try-error' %in% class(model)) b <- rep(NA, ", 
                   L, ") else { res=summary(model)$coefficients; b <- res[2,", colind.text, 
                   "]};invisible(b)}")
  return(script)
}

## Internal function: create iterator for bulk matrix by column

iblkcol_lag <- function(M, ...) {
  i <- 1
  it <- iterators::idiv(ncol(M), ...)
  
  nextEl <- function() {
    n <- iterators::nextElem(it)
    r <- seq(i, length = n)
    i <<- i + n
    M[, r, drop = FALSE]
  }
  obj <- list(nextElem = nextEl)
  class(obj) <- c("abstractiter", "iter")
  obj
}

## Internal function: scale data

scaleto <- function(dat) {
  if (is.null(dat)) 
    return(list(dn = NULL, d = NULL, ds = NULL))
  dat_scale <- scale(dat)
  dat_names <- names(dat)
  if (any(class(dat) %in% c("matrix", "data.frame", "data.table"))) {
    dat_names <- colnames(dat)
    dat <- as.matrix(data.frame(dat_scale))
  } else {
    dat_names <- names(dat)
    dat <- as.numeric(dat_scale)
  }
  dat_scale <- as.numeric(attributes(dat_scale)[["scaled:scale"]])
  return(list(dn = dat_names, d = dat, ds = dat_scale))
}

# Internal function: Sure Independent Screening

himasis <- function(Y, M, X, COV, glm.family, modelstatement, 
                    parallel = TRUE, ncore = 2, verbose = TRUE) {
  L.M <- ncol(M)
  M.names <- colnames(M)
  
  # if (scale) 
    # M_scale <- scaleto(M) else M_scale <- M
  
  X <- data.frame(X)
  X <- data.frame(model.matrix(~., X))[, -1]
  
  
  if (is.null(COV)) {
    datarun <- data.frame(Y = Y, Mone = NA, X = X)
    modelstatement <- modelstatement
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1]
    conf.names <- colnames(COV)
    datarun <- data.frame(Y = Y, Mone = NA, X = X, COV = COV)
    modelstatement <- eval(parse(text = (paste0(modelstatement, "+", 
                                                paste0(paste0("COV.", conf.names), collapse = "+")))))
  }
  
  doOne <- eval(parse(text = doOneGen(paste0("try(glm(modelstatement, family = ", 
                                             glm.family, ", data = datarun))"), "c(1,4)")))
  
  checkParallel("Sure Independent Screening", parallel, ncore, verbose)
  
  results <- foreach(n = iterators::idiv(L.M, chunks = ncore), 
                     M_chunk = iblkcol_lag(M, chunks = ncore), 
                     .combine = "cbind") %dopar% {sapply(seq_len(n), doOne, datarun, M_chunk)}

  colnames(results) <- M.names
  return(results)
}