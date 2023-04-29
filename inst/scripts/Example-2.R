# Generate HIMA example based on real-world dataset - 2 (logistic outcome)

set.seed(1029)
p <- 300 # the dimension of mediators
q <- 2    # the dimension of covariates
n <- 300  # sample size
alpha <- matrix(0,1,p) # the coefficients for X -> M
beta <- matrix(0,1,p) # the coefficients for M -> Y
alpha[1:5] <- 0.7
beta[1:5] <- 0.5
zeta <- matrix(0.01,p,q) # the coefficients of covariates for X -> M
eta <- matrix(0.01,1,q) # the coefficients of covariates for M -> Y
gamma <- 0.5 # the direct effect
gamma_total <- gamma + alpha%*%t(beta) # the total effect
X <- matrix(rbinom(n, size=1, prob=0.6),n,1) # expoure: 1=treatment; 0=placebo
Z <- matrix(0,n,2) # covariates
Z[,1] <- rbinom(n, size=1, prob=0.5) # sex:  male=1; female=0
Z[,2] <- sample(18:65,n,replace=TRUE)   # age ranging from 18-65
#Z[,2] <- (Z[,2]-mean(Z[,2]))/sd(Z[,2]) # scaled age, so we add a note that age is scaled before modeling
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

Example2 <- list(PhenoData = pheno, Mediator = M)

usethis::use_data(Example2, overwrite = TRUE)
