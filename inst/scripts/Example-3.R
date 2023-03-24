# Generate HIMA example based on real-world dataset - 3 (survival outcome)

set.seed(1355)
n <- 200 # sample size
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
beta[1]  <- 0.5
beta[2] <-  0.5
beta[3]  <- 0.5
alpha <- matrix(0,1,p)
alpha[1]  <- 0.5
alpha[2] <-  0.5
alpha[3]  <- 0.5
eta <- matrix(0.5,1,q)
zeta <- matrix(0.3,p,q)
gamma <- matrix(0.5,1,1)
X <- matrix(rbinom(n, size=1, prob=0.6),n,1) # expoure: 1=treatment; 0=placebo
Z <- matrix(0,n,2) # covariates
Z[,1] <- rbinom(n, size=1, prob=0.5) # sex:  male=1; female=0
Z[,2] <- sample(18:65,n,replace=TRUE)   # age ranging from 18-65
Z[,2] <- (Z[,2]-mean(Z[,2]))/sd(Z[,2]) # scaled age, so we add a note that age is scaled before modeling
mu <- matrix(0,p,1)
e <- mvrnorm(n, mu, sigma_e)  # the error terms
M <- matrix(0,n,p)
M <- X%*%(alpha) + Z%*%t(zeta) + e # the mediator matrix
MZ <- cbind(M,Z,X)
beta_gamma <- cbind(beta,eta,gamma)
u <- runif(n, 0, 1)
T <- matrix(0,n,1)
for (i in 1:n){
  T[i] <- -log(1-u[i])*exp(-sum(beta_gamma*MZ[i,]))
}
C <- runif(n, min = 0, max = 50)  # generate censoring time
status <-  T < C  # the censoring status
censoring_rate  <-  1 -  mean(as.numeric(status))  #  censoring rate is about 20%
OT <- apply(cbind(C,T),1,min) # the observed failure times

colnames(M) <- paste0("M", 1:100)
rownames(M) <- paste0("S", 1:200)

pheno <- data.frame(Treatment = X, Disease = status, Time = OT, Sex = Z[,1], Age = Z[,2])

Example3 <- list(PhenoData = pheno, Mediator = M)

usethis::use_data(Example3, overwrite = TRUE)
