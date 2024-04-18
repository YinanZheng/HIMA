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
              format(Sys.time(), "%X"), ")")
      doParallel::registerDoParallel(ncore)
  } else {
    if (verbose) 
      message("    Running ", program.name, " with single core...   (", 
              format(Sys.time(), "%X"), ")")
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



## Internal function: check data and scale 

dat <- data.frame(y = 1:12, x = c(rep("A",3), rep("B", 3), rep("C", 3), rep("D", 3)))
dat
f <- y~x

formula <- f
model.matrix(f, dat)

model.matrix(update(f, .~.-1), dat)


checkscale <- function(formula, dat) {
  if (is.null(dat)) 
    return(list(dn = NULL, d = NULL, ds = NULL))
  if(sum(is.na(dat))>0)
    return("Missing")
  
  
  
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
  
  if(glm.family == "gaussian")
  {
    doOne <- eval(parse(text = doOneGen(paste0("try(glm(modelstatement, family = ", 
                                               glm.family, ", data = datarun))"), "c(1,4)")))
  } else if(glm.family == "negbin") {
    doOne <- eval(parse(text = doOneGen(paste0("try(MASS::glm.nb(modelstatement, data = datarun))"), "c(1,4)")))
  } else {
    stop(paste0("Screening family ", glm.family, " is not supported."))
  }

  
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



# Internal function: null_estimation
# A function to estimate the proportions of the three component nulls
# This is from HDMT package version < 1.0.4

null_estimation <- function(input_pvalues, lambda = 0.5) {
  ## input_pvalues is a matrix with 2 columns of p-values, the first column is p-value for exposure-mediator association, the second column is p-value for mediator-outcome association adjusted for exposure
  ## lambda is the threshold for pi_{00} estimation, default 0.5
  #check input
  if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) !=2)
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues),nrow=nrow(input_pvalues))
  if (sum(stats::complete.cases(input_pvalues))<nrow(input_pvalues))
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[stats::complete.cases(input_pvalues),]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues)<1)
    stop("input_pvalues doesn't have valid p-values")
  
  pcut <- seq(0.1,0.8,0.1) 
  frac1 <- rep(0,8)
  frac2 <- rep(0,8)
  frac12<- rep(0,8)
  for (i in 1:8) {
    frac1[i] <- mean(input_pvalues[,1]>=pcut[i])/(1-pcut[i])
    frac2[i] <- mean(input_pvalues[,2]>=pcut[i])/(1-pcut[i]) 
    frac12[i]<- mean(input_pvalues[,2]>=pcut[i] & input_pvalues[,1]>=pcut[i])/(1-pcut[i])^2
  }  
  
  ## use the median estimates for pi00 ##
  
  alpha00 <- min(frac12[pcut==lambda],1)
  
  ## alpha1 is the proportion of nulls for first p-value 
  ## alpha2 is the proportion of nulls for second p-value 
  
  if (stats::ks.test(input_pvalues[,1],"punif",0,1,alternative="greater")$p>0.05) alpha1 <- 1 else   alpha1 <- min(frac1[pcut==lambda],1)  
  if (stats::ks.test(input_pvalues[,2],"punif",0,1,alternative="greater")$p>0.05) alpha2 <- 1 else   alpha2 <- min(frac2[pcut==lambda],1)
  
  
  if (alpha00==1) {
    alpha01 <- 0
    alpha10 <- 0
    alpha11 <- 0
  } else {    
    if (alpha1==1  & alpha2==1) {
      alpha01 <- 0
      alpha10 <- 0
      alpha11 <- 0
      alpha00 <- 1
    }  
    
    if (alpha1==1  & alpha2!=1) {
      alpha10 <- 0
      alpha11 <- 0
      alpha01 <- alpha1-alpha00
      alpha01 <- max(0,alpha01)
      alpha00 <- 1-alpha01
    }  
    
    if (alpha1!=1  & alpha2==1) {
      alpha01 <- 0
      alpha11 <- 0
      alpha10 <- alpha2-alpha00
      alpha10 <- max(0,alpha10)
      alpha00 <- 1-alpha10
    }  
    
    if (alpha1!=1  & alpha2!=1) {
      alpha10 <- alpha2-alpha00
      alpha10 <- max(0,alpha10)
      alpha01 <- alpha1-alpha00
      alpha01 <- max(0,alpha01)
      
      if ((1-alpha00-alpha01-alpha10)<0) {
        alpha11 <- 0
        alpha10 <- 1- alpha1
        alpha01 <- 1- alpha2
        alpha00 <- 1- alpha10 - alpha01
      }  else {
        alpha11 <-  1-alpha00-alpha01-alpha10
      }  
    }  
  }
  alpha.null <- list(alpha10=alpha10,alpha01=alpha01,alpha00=alpha00,alpha1=alpha1,alpha2=alpha2)
  return(alpha.null)
}



# Internal function: DLASSO_fun
# A function perform de-biased lasso estimator used by function "microHIMA"

DLASSO_fun <- function(X,Y){
  n <- dim(X)[1]
  p <- dim(X)[2]
  fit = glmnet(X, Y, alpha = 1)
  cv.fit <- cv.glmnet(X, Y, alpha = 1)
  beta_0 <- coef(fit, s = cv.fit$lambda.min)[2:(p+1)]
  #
  fit <- glmnet(X[,2:p], X[,1], alpha = 1)
  cv.fit <- cv.glmnet(X[,2:p], X[,1], alpha = 1)
  phi_hat <- coef(fit, s = cv.fit$lambda.min)[2:p]
  ##
  R <- X[,1] - X[,2:p]%*%t(t(phi_hat))
  E <- Y - X%*%t(t(beta_0))
  beta_1_hat <- beta_0[1] + sum(R*E)/sum(R*X[,1]) #  The de-biased lasso estimator
  ##
  sigma_e2 <- sum(E^2)/(n-length(which(beta_0!=0)))
  
  sigma_beta1_hat <- sqrt(sigma_e2)*sqrt(sum(R^2))/abs(sum(R*X[,1]))
  
  results <- c(beta_1_hat,sigma_beta1_hat)
  return(results)
}



# Internal function: rdirichlet
# A function generate random number from Dirichlet distribution.

rdirichlet <- function (n = 1, alpha) 
{
  Gam <- matrix(0, n, length(alpha))
  for (i in 1:length(alpha)) Gam[, i] <- stats::rgamma(n, shape = alpha[i])
  Gam/rowSums(Gam)
}



# Internal function: DACT
# A function to perform Divide-Aggregate Composite-null Test (DACT) by Liu et al. (2020).
# p value is corrected by JC method Jin and Cai (2007).
# This function is used in eHIMA

DACT = function(p_a,p_b){
  Z_a = stats::qnorm(p_a,lower.tail = F)
  Z_b = stats::qnorm(p_b,lower.tail = F)
  pi0a = 1 - .nonnullPropEst(Z_a,0,1)
  pi0b = 1 - .nonnullPropEst(Z_b,0,1)
  if(pi0a > 1){
    pi0a = 1
  }
  if(pi0b >1){
    pi0b = 1
  }
  p.mat = cbind(p_a,p_b)
  p3 = (apply(p.mat,1,max))^2
  wg1 = pi0a*(1-pi0b)
  wg2 = (1-pi0a)*pi0b
  wg3 = pi0a*pi0b
  wg.sum = wg1 + wg2 + wg3
  wg.std = c(wg1,wg2,wg3)/wg.sum
  p_dact = wg.std[1]*p_a + wg.std[2]*p_b + wg.std[3]*p3
  p_dact = .JCCorrect(p_dact)
  return(p_dact)
}

.JCCorrect = function(pval){
  z = stats::qnorm(pval,lower.tail = F)
  res= .nullParaEst(z)
  pval.JC = stats::pnorm(z,mean = res$mu,sd = res$s,lower.tail = F)
  return(pval.JC)
}

.nonnullPropEst <- function(x,u,sigma)
{
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)
  
  epsest=NULL
  
  for (j in 1:length(tt)) {
    
    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi
    
    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    }
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}

.nullParaEst<-function (x,gamma=0.1)
{
  n = length(x)
  t = c(1:1000)/200
  
  gan    = n^(-gamma)
  that   = 0
  shat   = 0
  uhat   = 0
  epshat = 0
  
  phiplus   = rep(1,1000)
  phiminus  = rep(1,1000)
  dphiplus  = rep(1,1000)
  dphiminus = rep(1,1000)
  phi       = rep(1,1000)
  dphi      = rep(1,1000)
  
  for (i in 1:1000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
  }
  
  ind = min(c(1:1000)[(phi - gan) <= 0])
  tt = t[ind]
  a  = phiplus[ind]
  b  = phiminus[ind]
  da = dphiplus[ind]
  db = dphiminus[ind]
  c  = phi[ind]
  
  that   = tt
  shat   = -(a*da + b*db)/(tt*c*c)
  shat   = sqrt(shat)
  uhat   = -(da*b - db*a)/(c*c)
  epshat = 1 - c*exp((tt*shat)^2/2)
  
  return(musigma=list(mu=uhat,s=shat))
}
