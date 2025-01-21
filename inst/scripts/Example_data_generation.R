# Generate HIMA example dataset based on real-world dataset

set.seed(881029)

######################################################
# Dataset - 1 (linear outcome)
######################################################

p <- 300 # the dimension of mediators
q <- 2    # the dimension of covariates
n <- 300  # sample size
alpha <- matrix(0,1,p) # the coefficients for X -> M
beta <- matrix(0,1,p) # the coefficients for M -> Y
alpha[1:3] <- c(0.7, 0.8, -0.9)
beta[1:3] <- c(0.9, 0.8, -0.7)
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

p <- 300 # the dimension of mediators
q <- 2    # the dimension of covariates
n <- 300  # sample size
alpha <- matrix(0,1,p) # the coefficients for X -> M
beta <- matrix(0,1,p) # the coefficients for M -> Y
alpha[1:3] <- c(0.7, 0.8, -0.9)
beta[1:3] <- c(0.9, 0.8, -0.7)
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

n <- 200 # sample size
p <- 200 # the dimension of mediators
q <- 2    # the dimenson of covariates
sigma_e <- matrix(0,p,p)
for (i in 1:p){
  for (j in 1:p){
    sigma_e[i,j] <- 0.25^{abs(i-j)}
    sigma_e[j,i] <- sigma_e[i,j]
  }
}
beta <- matrix(0,1,p)
beta[1]  <- 0.7
beta[2] <-  0.8
beta[3]  <- -0.9
alpha <- matrix(0,1,p)
alpha[1]  <- 0.9
alpha[2] <-  0.8
alpha[3]  <- -0.7
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
C <- runif(n, min = 0, max = 1)  # generate censoring time
status <-  T < C  # the censoring status
censoring_rate  <-  1 -  mean(as.numeric(status))  #  censoring rate is about 20%
OT <- apply(cbind(C,T),1,min) # the observed failure times

colnames(M) <- paste0("M", 1:p)
rownames(M) <- paste0("S", 1:n)

pheno <- data.frame(Treatment = X, Status = status, Time = OT, Sex = Z[,1], Age = Z[,2])

SurvivalData <- list(PhenoData = pheno, Mediator = M)
usethis::use_data(SurvivalData, overwrite = TRUE)



######################################################
# Dataset - 4 (microbiome mediators)
######################################################

n <- 100 # the sample size
p <- 100 # the dimension of mediator
a <- matrix(0,1,p) # the coefficient a in Eq. (16)
a[1] <- 1/3
a[2] <- 1/4
a[3] <- 1/5
eta <- matrix(0.01,2,1) # the coefficients of covariates for M -> Y
UR <- runif(p-3,0,1) # for generating a_4,...., a_p
a[4:p] <- (1-sum(a[1:3]))*UR/sum(UR)
b <- matrix(0,1,p) # the coefficient b in Eq. (17)
b[1] <- 1.3
b[2] <- -0.7
b[3] <- -0.6
SIG <- diag(p-1) + matrix(1,p-1,1)%*%matrix(1,1,p-1) # the covariance of e
SIG_e <- 2*SIG
M0 <- matrix(runif(n*p, min = 0, max = 1),n,p) # the baseline composition m0
for (i in 1:n){
  M0[i,] <- M0[i,]/sum(M0[i,])
}
trt <- matrix(rbinom(n, size=1, prob=0.6),n,1) # trt: 1=treatment; 0=placebo
Z <- matrix(0,n,2) # covariates
Z[,1] <- rbinom(n, size=1, prob=0.5) # sex:  male=1; female=0
Z[,2] <- sample(18:65,n,replace=TRUE)   # age ranging from 18-65
X <- trt + Z%*%eta  #  note that X is not the exposure, but only for generating M
e <- MASS::mvrnorm(n, rep(0, p-1), SIG_e)
U1 <- matrix(0,n,p) # generating compositional error e in Eq. (16)
for (i in 1:n){
  U1[i, 1:(p-1)] <- exp(e[i,])/(1+sum(exp(e[i,])))
  U1[i, p]       <- 1/(1+sum(exp(e[i,])))
}
aX <- matrix(0,n,p) # the a^X term
for (i in 1:n){
  aX[i,] <- a^(X[i])/sum(a^(X[i]))
}
M_raw <- matrix(0,n,p) # the compositional mediators
for (i in 1:n){
  A1 <- M0[i,]*aX[i,]/sum(M0[i,]*aX[i,]) # the definition of  perturbation operator for M
  M_raw[i,] <- A1*U1[i,]/sum(A1*U1[i,]) # the definition of  perturbation operator for M
}
U2 <- rnorm(n,0,2) # the error term for Y in Eq. (17)
Y <- 1+ 0.5*trt + log(M_raw)%*%t(b) + Z%*%(eta)+ U2 # the outcome Y
OTU <- M_raw # high-dimensional compositional matrix

colnames(OTU) <- paste0("M", 1:p)
rownames(OTU) <- paste0("S", 1:n)

pheno <- data.frame(Treatment = trt, Outcome = Y, Sex = Z[,1], Age = Z[,2])

MicrobiomeData <- list(PhenoData = pheno, Mediator = OTU)
usethis::use_data(MicrobiomeData, overwrite = TRUE)



######################################################
# Dataset - 5 (quantile mediation analysis)
######################################################

p <- 300 # the dimension of mediators
q <- 2    # the dimension of covariates
n <- 300  # sample size
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
