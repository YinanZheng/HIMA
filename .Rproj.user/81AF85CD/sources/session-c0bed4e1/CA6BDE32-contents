# Generate HIMA example dataset based on real-world dataset



######################################################
# Dataset - 1 (linear outcome)
######################################################
set.seed(1837)
p <- 100 # the dimension of mediators
q <- 2    # the dimension of covariates
n <- 500  # sample size
alpha <- matrix(0,1,p) # the coefficients for X -> M
beta <- matrix(0,1,p) # the coefficients for M -> Y
alpha[1:3] <- c(0.5, 0.8, -1)
beta[1:3] <- c(1, 0.8, -0.8)
zeta <- matrix(0.01,p,q) # the coefficients of covariates for X -> M
eta <- matrix(0.01,1,q) # the coefficients of covariates for M -> Y
gamma <- 1 # the direct effect
gamma_total <- gamma + alpha%*%t(beta) # the total effect
X <- matrix(rbinom(n, size=1, prob=0.6),n,1) # expoure: 1=treatment; 0=placebo
Z <- matrix(0,n,2) # covariates
Z[,1] <- rbinom(n, size=1, prob=0.5) # sex:  male=1; female=0
Z[,2] <- sample(18:65,n,replace=TRUE)   # age ranging from 18-65
e <- matrix(rnorm(n*p, mean = 0, sd = 1.5),n,p)
E <- matrix(rnorm(n, mean = 0, sd = 1),n,1)
M <-  X%*%(alpha) + Z%*%t(zeta) + e # the mediator matrix
Y <-  X*gamma + M%*%t(beta) + Z%*%t(eta) + E # the response Y

colnames(M) <- paste0("M", 1:p)
rownames(M) <- paste0("S", 1:n)

pheno <- data.frame(Treatment = X, Outcome = Y, Sex = Z[,1], Age = Z[,2])

ContinuousOutcome <- list(PhenoData = pheno, Mediator = M)
usethis::use_data(ContinuousOutcome, overwrite = TRUE)





######################################################
# Dataset - 2 (logistic outcome)
######################################################
set.seed(185)
p <- 100 # the dimension of mediators
q <- 2    # the dimension of covariates
n <- 500  # sample size
alpha <- matrix(0,1,p) # the coefficients for X -> M
beta <- matrix(0,1,p) # the coefficients for M -> Y
alpha[1:3] <- c(0.8, 0.8, -0.9)
beta[1:3] <- c(0.7, 0.7, -0.9)
zeta <- matrix(0.01,p,q) # the coefficients of covariates for X -> M
eta <- matrix(0.01,1,q) # the coefficients of covariates for M -> Y
gamma <- 1 # the direct effect
gamma_total <- gamma + alpha%*%t(beta) # the total effect
X <- matrix(rbinom(n, size=1, prob=0.6),n,1) # expoure: 1=treatment; 0=placebo
Z <- matrix(0,n,2) # covariates
Z[,1] <- rbinom(n, size=1, prob=0.5) # sex:  male=1; female=0
Z[,2] <- sample(18:65,n,replace=TRUE)   # age ranging from 18-65
e <- matrix(rnorm(n*p, mean = 0, sd = 1.5),n,p)
M <-  X%*%(alpha) + Z%*%t(zeta) + e # the mediators
pi <-   1 - 1/(1 + exp(X*gamma + M%*%t(beta) + Z%*%t(eta)))
Y <- matrix(0,n,1) # the binary Y
for (j in 1:n){
  Y[j] <- rbinom(1, 1, pi[j])
}

colnames(M) <- paste0("M", 1:p)
rownames(M) <- paste0("S", 1:n)

pheno <- data.frame(Treatment = X, Disease = Y, Sex = Z[,1], Age = Z[,2])

BinaryOutcome <- list(PhenoData = pheno, Mediator = M)
usethis::use_data(BinaryOutcome, overwrite = TRUE)





######################################################
# Dataset - 3 (survival outcome)
######################################################
set.seed(18675)
n <- 300 # sample size
p <- 100 # the dimension of mediators
q <- 2    # the dimenson of covariates
sigma_e <- matrix(0,p,p)
for (i in 1:p){
  for (j in 1:p){
    sigma_e[i,j] <- 0.25^{abs(i-j)}
    sigma_e[j,i] <- sigma_e[i,j]
  }
}
beta <- matrix(0,1,p)
beta[1]  <- 0.8
beta[2] <-  0.8
beta[3]  <- -0.9
alpha <- matrix(0,1,p)
alpha[1]  <- 0.9
alpha[2] <-  0.9
alpha[3]  <- -0.8
eta <- matrix(0.01,1,q)
zeta <- matrix(0.01,p,q)
gamma <- matrix(1,1,1)
X <- matrix(rbinom(n, size=1, prob=0.6),n,1) # expoure: 1=treatment; 0=placebo
Z <- matrix(0,n,2) # covariates
Z[,1] <- rbinom(n, size=1, prob=0.5) # sex:  male=1; female=0
Z[,2] <- sample(18:65,n,replace=TRUE)   # age ranging from 18-65
mu <- matrix(0,p,1)
e <- MASS::mvrnorm(n, mu, sigma_e)  # the error terms
M <- matrix(0,n,p)
M <- X%*%(alpha) + Z%*%t(zeta) + e # the mediator matrix
MZ <- cbind(M,Z,X)
beta_gamma <- cbind(beta,eta,gamma)
u <- runif(n, 0, 1)
T <- matrix(0,n,1)
for (i in 1:n){
  T[i] <- -log(1-u[i])*exp(-sum(beta_gamma*MZ[i,]))
}
C <- runif(n, min = 0, max = 1.5)  # generate censoring time
status <-  T < C  # the censoring status
censoring_rate  <-  1 -  mean(as.numeric(status))  #  censoring rate is about 20%
OT <- apply(cbind(C,T),1,min) # the observed failure times

colnames(M) <- paste0("M", 1:p)
rownames(M) <- paste0("S", 1:n)

pheno <- data.frame(Treatment = X, Status = status, Time = OT, Sex = Z[,1], Age = Z[,2])

SurvivalData <- list(PhenoData = pheno, Mediator = M)
usethis::use_data(SurvivalData, overwrite = TRUE)





######################################################
# Dataset - 4 (longitudinal mediator + survival outcome)
######################################################
set.seed(18675)

n <- 300   # sample size
p <- 100   # the dimension of mediators
q <- 2     # the dimension of covariates

alpha <- c(c(0.5, 0.4, 0.3),
           c(0.5, 0.4, 0.3),
           rep(0, 3),
           rep(0, p - 9))
eta    <- matrix(0.01, q, p)
theta1 <- 0.4
theta2 <- rep(0.01, q)
beta   <- c(c(0.5, 0.4, 0.3),
            rep(0, 3),
            c(0.5, 0.4, 0.3),
            rep(0, p - 9))

a1 <- rnorm(n, 0, 0.2)

## ---------------------------------------------------
## 1. Generate Sex / Age, then standardize for modeling
## ---------------------------------------------------

# Raw covariates (to be saved in PhenoData)
Sex_raw <- rbinom(n, size = 1, prob = 0.5)         # sex: 1 = male, 0 = female
Age_raw <- sample(18:65, n, replace = TRUE)        # age 18–65

# Standardized version used in the longitudinal model
Z <- scale(cbind(Sex_raw, Age_raw))                # mean 0, sd 1
colnames(Z) <- c("Sex_z", "Age_z")

## ---------------------------------------------------
## 2. Exposure and time setup
## ---------------------------------------------------

x <- rnorm(n, 0, 0.5)

# intervals D1 = [0, t1), D2 = [t1, +infty)
interval.max <- 2
delta1 <- 0.2
t <- 0:(interval.max - 1) * delta1

beta0.ti <- c(theta1, theta2)
beta0.tv <- beta
beta0    <- c(beta0.ti, rep(beta0.tv, each = interval.max))
p.ti     <- 1 + q
p.tv     <- p
p.x      <- p.ti + p.tv * interval.max

## ---------------------------------------------------
## 3. Longitudinal mediators
## ---------------------------------------------------

M <- matrix(0, n, p * interval.max)
for (k in 1:p) {
  b1 <- rnorm(n, 0, 0.2)
  for (j in 1:interval.max) {
    M[, (k - 1) * interval.max + j] <-
      as.vector(x * alpha[k] + Z %*% eta[, k]) +
      a1 + b1 + rnorm(n, 0, 1)
  }
}

########## Austin Weibull #################
delta2 <- 1
scale  <- 1
nu     <- 4

# Design matrix used for hazard
xx0 <- cbind(x, Z, M, 1:n)   # last column is subject id
u   <- runif(n, 0, 1)

H    <- matrix(0, nrow = n, ncol = interval.max)
R    <- matrix(0, nrow = n, ncol = interval.max)
Hinv <- matrix(0, nrow = n, ncol = interval.max)
flag <- matrix(0, nrow = n, ncol = interval.max)

for (i in 1:interval.max) {
  ind <- c(1:p.ti,
           p.ti + (0:(p.tv - 1)) * interval.max + i)
  H[, i] <- scale * exp(xx0[, ind, drop = FALSE] %*%
                          c(beta0.ti, beta0.tv))
  if (i > 1) {
    R[, i]    <- R[, i - 1] + H[, i - 1] * (t[i]^nu - t[i - 1]^nu)
    flag[, i] <- -log(u) < R[, i]
  }
  Hinv[, i] <- ((-log(u) - R[, i]) / H[, i] + t[i]^nu)^(1 / nu)
  Hinv[, i][is.nan(Hinv[, i])] <- Inf
}

h.ind <- interval.max - rowSums(flag)

## ---------------------------------------------------
## 4. Build counting-process style data
## ---------------------------------------------------

xx <- NULL

for (i in 1:interval.max) {
  
  ind.col  <- c(1:p.ti,
                p.ti + (0:(p.tv - 1)) * interval.max + i,
                ncol(xx0))
  ind.elig <- h.ind >= i
  n.elig   <- sum(ind.elig)
  
  if (n.elig == 0) next
  
  t0 <- rep(t[i], n.elig)
  t1 <- Hinv[ind.elig, i]
  
  if (i == 1) {
    t1[t1 > t[i] + delta1] <- t[i] + delta1
    status <- t1 == Hinv[ind.elig, i]
    
    xx <- cbind(
      t0     = t0,
      t1     = t1,
      status = status,
      xx0[ind.elig, ind.col, drop = FALSE]
    )
  } else {
    t1[t1 > t[i] + delta2] <- t[i] + delta2
    status <- t1 == Hinv[ind.elig, i]
    
    xx <- rbind(
      xx,
      cbind(
        t0     = t0,
        t1     = t1,
        status = status,
        xx0[ind.elig, ind.col, drop = FALSE]
      )
    )
  }
}

data <- xx[order(xx[, ncol(xx)]), ]

id     <- data[, ncol(data)]
tstart <- data[, "t0"]
tstop  <- data[, "t1"]
status <- data[, "status"]
X      <- data[, 4]

## ---------------------------------------------------
## 5. PhenoData: attach raw Sex / Age by id
## ---------------------------------------------------

PhenoData <- data.frame(
  ID       = id,
  Tstart   = tstart,
  Tstop    = tstop,
  Status   = status,
  Treatment = X,
  Sex      = Sex_raw[id],
  Age      = Age_raw[id],
  row.names = NULL
)

## ---------------------------------------------------
## 6. Mediator matrix (time-varying part)
## ---------------------------------------------------

Mt <- data[, (4 + q + 1):(4 + q + p), drop = FALSE]
Mediator <- as.matrix(Mt)
colnames(Mediator) <- paste0("M", 1:ncol(Mediator))

SurvivalLongData <- list(
  PhenoData = PhenoData,
  Mediator  = Mediator
)

usethis::use_data(SurvivalLongData, overwrite = TRUE)





######################################################
# Dataset - 6 (quantile mediation analysis)
######################################################
set.seed(1753)
p <- 100 # the dimension of mediators
q <- 2    # the dimension of covariates
n <- 500  # sample size
alpha <- matrix(0,1,p) # the coefficients for X -> M
beta <- matrix(0,1,p) # the coefficients for M -> Y

alpha[1] <- 0.7
alpha[2] <- 0.8
alpha[3] <- -0.9
beta[1] <- 0.9
beta[2] <- 0.8
beta[3] <- -0.7

sigma_e <- matrix(0,p,p)
rou <- 0.25  # the correlation of  M
for (i in 1:p) {
  for (j in 1:p) {
    sigma_e[i,j]=(rou^(abs(i-j)));
  }
}

E <- matrix(rt(n, 3)) # E is t distribution

zeta <- matrix(0.01,p,q) # the coefficients of covariates for X -> M
eta <- matrix(0.01,1,q) # the coefficients of covariates for M -> Y
gamma <- 0.8 # the direct effect
gamma_total <- gamma + alpha%*%t(beta) # the total effect
X <- matrix(rnorm(n, mean = 0, sd = 2),n,1) # expoure: 1=treatment; 0=placebo
Z <- matrix(0,n,q) # covariates
Z[,1] <- rbinom(n, size=1, prob=0.5) # sex:  male=1; female=0
Z[,2] <- sample(18:65,n,replace=TRUE)   # age ranging from 18-65

mu <- matrix(0,p,1)
e <- MASS::mvrnorm(n, mu, sigma_e)
M <-  X%*%(alpha) + Z%*%t(zeta) + e # the mediators
Y <-  X*gamma + M%*%t(beta) + Z%*%t(eta) + E # the response Y

colnames(M) <- paste0("M", 1:p)
rownames(M) <- paste0("S", 1:n)

pheno <- data.frame(Treatment = X, Outcome = Y, Sex = Z[,1], Age = Z[,2])

QuantileData <- list(PhenoData = pheno, Mediator = M)
usethis::use_data(QuantileData, overwrite = TRUE)
