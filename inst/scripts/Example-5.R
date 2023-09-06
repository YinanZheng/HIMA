# Generate HIMA example based on real-world dataset - 5 (quantile mediation analysis)

library(quantreg)
library(conquer)
library(MASS)

p <- 500 # the dimension of mediators
q <- 2    # the dimension of covariates
n <- 300  # sample size
alpha <- matrix(0,1,p) # the coefficients for X -> M
beta <- matrix(0,1,p) # the coefficients for M -> Y

alpha[1] <- 0.7
alpha[3] <- 0.8
alpha[5] <- 0.9
beta[1] <- 0.9
beta[3] <- 0.8
beta[5] <- 0.7

sigma_e <- matrix(0,p,p)
rou <- 0.25  # the correlation of  M
for (i in 1:p) {
  for (j in 1:p) {
    sigma_e[i,j]=(rou^(abs(i-j)));
  }
}

E <- matrix(rt(n, 3)) # E is t distribution

zeta <- matrix(0.3,p,q) # the coefficients of covariates for X -> M
eta <- matrix(0.3,1,q) # the coefficients of covariates for M -> Y
gamma <- 1 # the direct effect
gamma_total <- gamma + alpha%*%t(beta) # the total effect
X <- matrix(rbinom(n, size=1, prob=0.6),n,1) # expoure: 1=treatment; 0=placebo
Z <- matrix(0,n,q) # covariates
Z[,1] <- rbinom(n, size=1, prob=0.5) # sex:  male=1; female=0
Z[,2] <- sample(18:65,n,replace=TRUE)   # age ranging from 18-65

mu <- matrix(0,p,1)
e <- mvrnorm(n, mu, sigma_e)
M <-  X%*%(alpha) + Z%*%t(zeta) + e # the mediators
Y <-  X*gamma + M%*%t(beta) + Z%*%t(eta) + E # the response Y

colnames(M) <- paste0("M", 1:p)
rownames(M) <- paste0("S", 1:n)

pheno <- data.frame(Treatment = X, Outcome = Y, Sex = Z[,1], Age = Z[,2])

Example5 <- list(PhenoData = pheno, Mediator = M)

usethis::use_data(Example5, overwrite = TRUE)
