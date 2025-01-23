# Internal function: parallel computing check

checkParallel <- function(program.name, parallel, ncore, verbose) {
  if (parallel && (ncore > 1)) {
    if (ncore > parallel::detectCores()) {
      message(
        "You requested ", ncore, " cores. There are only ",
        parallel::detectCores(), " in your machine!"
      )
      ncore <- parallel::detectCores()
    }
    if (verbose) {
      message(
        "    Running ", program.name, " with ", ncore, " cores in parallel...   (",
        format(Sys.time(), "%X"), ")"
      )
    }
    doParallel::registerDoParallel(ncore)
  } else {
    if (verbose) {
      message(
        "    Running ", program.name, " with single core...   (",
        format(Sys.time(), "%X"), ")"
      )
    }
    registerDoSEQ()
  }
}



## Internal function: doOne code generater

doOneGen <- function(model.text, colind.text) {
  L <- length(eval(parse(text = colind.text)))
  script <- paste0(
    "doOne <- function(i, datarun, Ydat){datarun$Mone <- Ydat[,i]; model <- ",
    model.text, ";if('try-error' %in% class(model)) b <- rep(NA, ",
    L, ") else { res=summary(model)$coefficients; b <- res[2,", colind.text,
    "]};invisible(b)}"
  )
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



# Internal function: Sure Independent Screening for hima

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
    COV <- data.frame(model.matrix(~., COV))[, -1, drop = FALSE]
    conf.names <- colnames(COV)
    if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    datarun <- data.frame(Y = Y, Mone = NA, X = X, COV)
    modelstatement <- eval(parse(text = (paste0(
      modelstatement, "+",
      paste0(colnames(COV), collapse = "+")
    ))))
  }

  if (glm.family == "gaussian") {
    doOne <- eval(parse(text = doOneGen(paste0(
      "try(glm(modelstatement, family = ",
      glm.family, ", data = datarun))"
    ), "c(1,4)")))
  } else if (glm.family == "negbin") {
    doOne <- eval(parse(text = doOneGen(paste0("try(MASS::glm.nb(modelstatement, data = datarun))"), "c(1,4)")))
  } else {
    stop(paste0("Screening family ", glm.family, " is not supported."))
  }


  checkParallel(tag, parallel, ncore, verbose)

  results <- foreach(
    n = iterators::idiv(L.M, chunks = ncore),
    M_chunk = iblkcol_lag(M, chunks = ncore),
    .combine = "cbind"
  ) %dopar% {
    sapply(seq_len(n), doOne, datarun, M_chunk)
  }

  colnames(results) <- M.names
  return(results)
}



# Internal function: process_var
# Helper function to process variables

process_var <- function(var, scale) {
  if (!is.null(var)) {
    if (scale) {
      return(scale(var))
    } else {
      return(as.matrix(var))
    }
  } else {
    return(NULL)
  }
}



# Internal function: LDPE (optimized)
# the code of Liu han's JRSSB paper for high-dimensional Cox model
# ID: the index of interested parameter
# X: the covariates matrix with n by p
# OT: the observed time = min(T,C)
# status: the censoring indicator I(T <= C)

LDPE_func <- function(ID, X, OT, status) {
  set.seed(1029)

  # Dimensions
  coi <- ID
  d <- ncol(X)
  n <- nrow(X)

  # Initialize penalty factor
  PF <- rep(1, d)
  PF[ID] <- 1

  # Semi-penalized initial estimator
  fit <- glmnet(X, survival::Surv(OT, status), family = "cox", alpha = 1, standardize = FALSE, penalty.factor = PF)
  cv.fit <- cv.glmnet(X, survival::Surv(OT, status), family = "cox", alpha = 1, standardize = FALSE, penalty.factor = PF)
  betas <- as.vector(coef(fit, s = cv.fit$lambda.min))

  # Precompute sorted times and order
  stime <- sort(OT)
  otime <- order(OT)

  # Preallocate storage
  la <- numeric(n)
  lb <- matrix(0, n, d - 1)
  Hs <- matrix(0, d, d)

  # Iterate through observations
  for (i in seq_len(n)) {
    if (status[otime[i]] == 1) {
      ind <- which(OT >= stime[i])
      X_ind <- X[ind, , drop = FALSE] 
      tmp_exp <- as.vector(exp(X_ind %*% betas)) 

      S0 <- sum(tmp_exp)
      S1 <- colSums(X_ind * tmp_exp)  
      S2 <- t(X_ind) %*% (X_ind * tmp_exp)

      la[i] <- -(X[otime[i], coi] - S1[coi] / S0)

      if (coi == 1) {
        lb[i, ] <- -(X[otime[i], -1] - S1[-1] / S0)
      } else if (coi == d) {
        lb[i, ] <- -(X[otime[i], -d] - S1[-d] / S0)
      } else {
        lb[i, ] <- -(X[otime[i], -c(coi)] - S1[-c(coi)] / S0)
      }

      Hs <- Hs + (S0 * S2 - tcrossprod(S1)) / (n * S0^2)
    }
  }

  # De-biased Lasso step
  fit_res <- glmnet(lb, la, alpha = 1, standardize = FALSE, intercept = FALSE, lambda = sqrt(log(d) / n))
  what <- as.vector(coef(fit_res)[-1])

  # Final estimate and variance
  if (coi == 1) {
    S <- betas[coi] - (mean(la) - crossprod(what, colMeans(lb))) / (Hs[coi, coi] - crossprod(what, Hs[-1, coi]))
    var <- Hs[coi, coi] - crossprod(what, Hs[-1, coi])
  } else if (coi == d) {
    S <- betas[coi] - (mean(la) - crossprod(what, colMeans(lb))) / (Hs[coi, coi] - crossprod(what, Hs[-d, coi]))
    var <- Hs[coi, coi] - crossprod(what, Hs[-d, coi])
  } else {
    S <- betas[coi] - (mean(la) - crossprod(what, colMeans(lb))) / (Hs[coi, coi] - crossprod(what, Hs[-c(coi), coi]))
    var <- Hs[coi, coi] - crossprod(what, Hs[-c(coi), coi])
  }

  beta_est <- S
  beta_SE <- sqrt(1 / (n * var))

  result <- c(beta_est, beta_SE)
  return(result)
}


# Internal function: null_estimation
# A function to estimate the proportions of the three component nulls
# This is from HDMT package version < 1.0.4 (optimized here)

null_estimation <- function(input_pvalues, lambda = 0.5) {
  ## Validate input
  if (is.null(ncol(input_pvalues)) || ncol(input_pvalues) != 2) {
    stop("`input_pvalues` must be a matrix or data frame with exactly 2 columns.")
  }
  input_pvalues <- as.matrix(input_pvalues)
  if (anyNA(input_pvalues)) {
    warning("`input_pvalues` contains NAs, which will be removed.")
    input_pvalues <- input_pvalues[complete.cases(input_pvalues), ]
  }
  if (nrow(input_pvalues) == 0) {
    stop("`input_pvalues` does not contain valid rows.")
  }
  
  ## Precompute threshold values
  pcut <- seq(0.1, 0.8, 0.1)
  one_minus_pcut <- 1 - pcut
  one_minus_pcut_sq <- one_minus_pcut^2
  
  ## Calculate fractions using vectorized operations
  frac1 <- colMeans(outer(input_pvalues[, 1], pcut, `>=`)) / one_minus_pcut
  frac2 <- colMeans(outer(input_pvalues[, 2], pcut, `>=`)) / one_minus_pcut
  frac12 <- colMeans(outer(input_pvalues[, 2], pcut, `>=`) & outer(input_pvalues[, 1], pcut, `>=`)) / one_minus_pcut_sq
  
  ## Estimate alpha00
  alpha00 <- min(frac12[pcut == lambda], 1)
  
  ## Estimate alpha1 and alpha2
  alpha1 <- if (stats::ks.test(input_pvalues[, 1], "punif", 0, 1, alternative = "greater")$p > 0.05) 1 else min(frac1[pcut == lambda], 1)
  alpha2 <- if (stats::ks.test(input_pvalues[, 2], "punif", 0, 1, alternative = "greater")$p > 0.05) 1 else min(frac2[pcut == lambda], 1)
  
  ## Estimate other proportions
  if (alpha00 == 1) {
    alpha01 <- alpha10 <- alpha11 <- 0
  } else {
    alpha01 <- alpha1 - alpha00
    alpha10 <- alpha2 - alpha00
    alpha01 <- max(0, alpha01)
    alpha10 <- max(0, alpha10)
    
    if (alpha1 == 1 && alpha2 == 1) {
      alpha00 <- 1
      alpha01 <- alpha10 <- alpha11 <- 0
    } else if ((1 - alpha00 - alpha01 - alpha10) < 0) {
      alpha11 <- 0
      alpha10 <- 1 - alpha1
      alpha01 <- 1 - alpha2
      alpha00 <- 1 - alpha10 - alpha01
    } else {
      alpha11 <- 1 - alpha00 - alpha01 - alpha10
    }
  }
  
  ## Return results
  list(alpha10 = alpha10, alpha01 = alpha01, alpha00 = alpha00, alpha1 = alpha1, alpha2 = alpha2)
}


# Internal function: DLASSO_fun
# A function perform de-biased lasso estimator used by function "hima_microbiome"

DLASSO_fun <- function(X, Y) {
  set.seed(1029)
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  fit <- glmnet(X, Y, alpha = 1)
  cv.fit <- cv.glmnet(X, Y, alpha = 1)
  beta_0 <- coef(fit, s = cv.fit$lambda.min)[2:(p + 1)]
  #
  fit <- glmnet(X[, 2:p], X[, 1], alpha = 1)
  cv.fit <- cv.glmnet(X[, 2:p], X[, 1], alpha = 1)
  phi_hat <- coef(fit, s = cv.fit$lambda.min)[2:p]
  ##
  R <- X[, 1] - X[, 2:p] %*% t(t(phi_hat))
  E <- Y - X %*% t(t(beta_0))
  beta_1_hat <- beta_0[1] + sum(R * E) / sum(R * X[, 1]) #  The de-biased lasso estimator
  ##
  sigma_e2 <- sum(E^2) / (n - length(which(beta_0 != 0)))
  
  sigma_beta1_hat <- sqrt(sigma_e2) * sqrt(sum(R^2)) / abs(sum(R * X[, 1]))
  
  results <- c(beta_1_hat, sigma_beta1_hat)
  return(results)
}

# Internal function: rdirichlet
# A function generate random number from Dirichlet distribution.

rdirichlet <- function(n = 1, alpha) {
  Gam <- matrix(0, n, length(alpha))
  for (i in seq_along(alpha)) Gam[, i] <- stats::rgamma(n, shape = alpha[i])
  Gam / rowSums(Gam)
}



# Internal function: DACT (optimized)
# A function to perform Divide-Aggregate Composite-null Test (DACT) by Liu et al. (2020).
# p value is corrected by JC method Jin and Cai (2007).
# This function is used in hima_efficient

DACT <- function(p_a, p_b) {
  Z_a <- stats::qnorm(p_a, lower.tail = FALSE)
  Z_b <- stats::qnorm(p_b, lower.tail = FALSE)
  pi0a <- 1 - .nonnullPropEst(Z_a, 0, 1)
  pi0b <- 1 - .nonnullPropEst(Z_b, 0, 1)
  pi0a <- min(pi0a, 1)
  pi0b <- min(pi0b, 1)
  
  p3 <- (pmax(p_a, p_b))^2
  wg1 <- pi0a * (1 - pi0b)
  wg2 <- (1 - pi0a) * pi0b
  wg3 <- pi0a * pi0b
  wg.sum <- wg1 + wg2 + wg3
  wg.std <- c(wg1, wg2, wg3) / wg.sum
  
  p_dact <- wg.std[1] * p_a + wg.std[2] * p_b + wg.std[3] * p3
  p_dact <- .JCCorrect(p_dact)
  return(p_dact)
}

.JCCorrect <- function(pval) {
  z <- stats::qnorm(pval, lower.tail = FALSE)
  res <- .nullParaEst(z)
  stats::pnorm(z, mean = res$mu, sd = res$s, lower.tail = FALSE)
}

.nonnullPropEst <- function(x, u, sigma) {
  z <- (x - u) / sigma
  xi <- seq(0, 1, length.out = 101)
  tmax <- sqrt(log(length(x)))
  tt <- seq(0, tmax, 0.1)
  
  weights <- 1 - abs(xi)  # Weights based on xi
  epsest <- numeric(length(tt))  # Preallocate results
  
  for (j in seq_along(tt)) {
    t <- tt[j]
    f <- exp((t * xi)^2 / 2)
    co <- rowMeans(cos(outer(t * xi, z, `*`)))  # Calculate cosine terms
    epsest[j] <- 1 - sum(weights * f * co) / sum(weights)
  }
  
  max(epsest)
}

.nullParaEst <- function(x, gamma = 0.1) {
  n <- length(x)
  t <- seq(0.005, 5, length.out = 1000)  # Define range for smoother spacing
  
  gan <- n^(-gamma)
  
  # Compute cos(s * x) and sin(s * x) using outer for vectorization
  cos_vals <- outer(t, x, FUN = function(t, x) cos(t * x))  # Matrix of size (length(t), length(x))
  sin_vals <- outer(t, x, FUN = function(t, x) sin(t * x))  # Matrix of size (length(t), length(x))
  
  # Compute phi and its derivatives
  phiplus <- rowMeans(cos_vals)  # Mean of each row
  phiminus <- rowMeans(sin_vals)  # Mean of each row
  dphiplus <- -rowMeans(sweep(sin_vals, 2, x, `*`))  # Broadcasting x across columns
  dphiminus <- rowMeans(sweep(cos_vals, 2, x, `*`))  # Broadcasting x across columns
  phi <- sqrt(phiplus^2 + phiminus^2)  # Magnitude of phiplus and phiminus
  
  # Find the first index where phi - gan <= 0
  ind <- which((phi - gan) <= 0)[1]
  
  if (is.na(ind)) {
    stop("Unable to find a suitable index where phi - gan <= 0.")
  }
  
  # Extract values at the identified index
  tt <- t[ind]
  a <- phiplus[ind]
  b <- phiminus[ind]
  da <- dphiplus[ind]
  db <- dphiminus[ind]
  c <- phi[ind]
  
  # Compute final estimates
  shat <- sqrt(-(a * da + b * db) / (tt * c^2))
  uhat <- -(da * b - db * a) / (c^2)
  epshat <- 1 - c * exp((tt * shat)^2 / 2)
  
  list(mu = uhat, s = shat)
}
