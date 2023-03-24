# Generate HIMA example based on real-world dataset - 1 (linear outcome)

set.seed(1025)
p <- 10000 # the dimension of mediators
q <- 2    # the dimension of covariates
n <- 500  # sample size
alpha <- matrix(0,1,p) # the coefficients for X -> M
beta <- matrix(0,1,p) # the coefficients for M -> Y
alpha[1:5] <- 0.5
beta[1:5] <- 0.5
zeta <- matrix(0.3,p,q) # the coefficients of covariates for X -> M
eta <- matrix(0.5,1,q) # the coefficients of covariates for M -> Y
gamma <- 0.5 # the direct effect
gamma_total <- gamma + alpha%*%t(beta) # the total effect
X <- matrix(rbinom(n, size=1, prob=0.6),n,1) # expoure: 1=treatment; 0=placebo
Z <- matrix(0,n,2) # covariates
Z[,1] <- rbinom(n, size=1, prob=0.5) # sex:  male=1; female=0
Z[,2] <- sample(18:65,n,replace=TRUE)   # age ranging from 18-65
Z[,2] <- (Z[,2]-mean(Z[,2]))/sd(Z[,2]) # scaled age, so we add a note that age is scaled before modeling
e <- matrix(rnorm(n*p, mean = 0, sd = 1.5),n,p)
E <- matrix(rnorm(n, mean = 0, sd = 1),n,1)
M <-  X%*%(alpha) + Z%*%t(zeta) + e # the mediator matrix
Y <-  X*gamma + M%*%t(beta) + Z%*%t(eta) + E # the response Y

colnames(M) <- paste0("M", 1:10000)
rownames(M) <- paste0("S", 1:500)

pheno <- data.frame(Treatment = X, Outcome = Y, Sex = Z[,1], Age = Z[,2])

Example1 <- list(PhenoData = pheno, Mediator = M)

usethis::use_data(Example1, overwrite = TRUE)
