# Internal function: parallel computing check

checkParallel <- function(program.name, parallel, ncore, verbose) {
  if (parallel & (ncore > 1)) {
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



## Internal function: scale data (obsolete function)

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
# Global variables:
globalVariables("n")
globalVariables("M_chunk")

himasis <- function(Y, M, X, COV, glm.family, modelstatement, 
                    parallel, ncore, verbose, tag) {
  L.M <- ncol(M)
  M.names <- colnames(M)
  
  X <- data.frame(X)
  X <- data.frame(model.matrix(~., X))[, -1]
  
  if (is.null(COV)) {
    if (verbose) message("    No covariate is adjusted")
    datarun <- data.frame(Y = Y, Mone = NA, X = X)
    modelstatement <- modelstatement
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1]
    conf.names <- colnames(COV)
    if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    datarun <- data.frame(Y = Y, Mone = NA, X = X, COV = COV)
    modelstatement <- eval(parse(text = (paste0(modelstatement, "+", 
                                                paste0(paste0("COV.", conf.names), collapse = "+")))))
  }
  
  doOne <- eval(parse(text = doOneGen(paste0("try(glm(modelstatement, family = ", 
                                             glm.family, ", data = datarun))"), "c(1,4)")))
  
  checkParallel(tag, parallel, ncore, verbose)
  
  results <- foreach(n = iterators::idiv(L.M, chunks = ncore), 
                     M_chunk = iblkcol_lag(M, chunks = ncore), 
                     .combine = "cbind") %dopar% {sapply(seq_len(n), doOne, datarun, M_chunk)}

  colnames(results) <- M.names
  return(results)
}



# Internal function: LDPE
# the code of Liu han's JRSSB paper for high-dimensional Cox model
# ID: the index of interested parameter
# X: the covariates matrix with n by p
# OT: the observed time = min(T,C)
# status: the censoring indicator I(T <= C)

LDPE_func <- function(ID, X, OT, status){
  coi <- ID
  x <- X
  d <- dim(x)[2]
  n <- dim(x)[1]

  ##Set of tuning parameters
  PF <- matrix(1,1,d)
  PF[ID] <- 1
  fit <- glmnet(x, survival::Surv(OT, status), family="cox", alpha = 1, standardize = FALSE,penalty.factor=PF)
  cv.fit <- cv.glmnet(x, survival::Surv(OT, status), family="cox", alpha = 1, standardize = FALSE,penalty.factor=PF)
  betas   <-   coef(fit, s = cv.fit$lambda.min)[1:d]  # the semi-penalized initial estimator  # initial estimator
  
  stime = sort(OT)          # Sorted survival/censored times
  otime = order(OT)         # Order of time
  
  Vs  = matrix(rep(0,d*d),nrow = d)
  Hs  = Vs                                 # Hessian
  ind = 0
  
  la  = rep(0,n)                           # Gradient w.r.t parameter of interest
  lb  = matrix(rep(0,(d-1)*n),nrow = n)    # Gradient w.r.t nuisance parameter (theta)
  i   = 1
  
  while( i<=n)
  {
    if (status[otime[i]]==1)
    {
      ind = which(OT >= stime[i])
      S0  = 0
      S1  = rep(0,d)
      S2  = matrix(rep(0,d*d),nrow = d)
      
      if (length(ind)>0)
      {
        for (j in 1:length(ind))
        {
          tmp = exp(x[ind[j],]%*%betas)
          S0  = S0 + tmp
          
          S1  = S1 + tmp %*%t(x[ind[j],])
          
          tmp = apply(tmp,1,as.numeric)     
          S2  = S2 + tmp*x[ind[j],]%*%t(x[ind[j],])
        }
      }
      S0 = apply(S0,1,as.numeric)
      
      la[i]  = -(x[otime[i],coi] - S1[coi]/S0)
      if (coi == 1)
      {
        lb[i,] = -(x[otime[i],c((coi+1):d)] - S1[c((coi+1):d)]/S0)
      } else if (coi == d){
        lb[i,] = -(x[otime[i],c(1:(coi-1))] - S1[c(1:(coi-1))]/S0)
      } else {
        lb[i,] = -(x[otime[i],c(1:(coi-1), (coi+1):d)] - S1[c(1:(coi-1), (coi+1):d)]/S0)
      }
      V   = S0*S2 - t(S1)%*%(S1)
      Hs  = Hs + V/(n*S0^2)          
    }
    i = i + 1
  }
  
  fit <- glmnet(lb,la,alpha = 1, standardize = FALSE,intercept = FALSE,lambda = sqrt(log(d)/n))
  what <- as.numeric(coef(fit)[2:d])   
  
  if (coi == 1)
  {
    S = betas[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c((coi+1):d),coi])
    var   = Hs[coi,coi] - t(what)%*%Hs[c((coi+1):d),coi]
  } else if (coi == d){
    S = betas[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1)),coi])
    var   = Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1)),coi]
  } else {
    S = betas[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1),(coi+1):d),coi])
    var   = Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1),(coi+1):d),coi]
  }
  
  beta_est <- S
  beta_SE  <- sqrt(1/(n*var))
  
  result <- c(beta_est,beta_SE)
  
  return(result)
}
