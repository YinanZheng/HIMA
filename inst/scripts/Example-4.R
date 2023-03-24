# Generate HIMA example based on real-world dataset - 4 (microbiome mediators)

set.seed(130)
n <- 200 # the sample size
p <- 50 # the dimension of mediator
a <- matrix(0,1,p) # the coefficient a in Eq. (16)
a[1] <- 1/3
a[2] <- 1/4
a[3] <- 1/5
eta <- matrix(0.5,2,1) # the coefficients of covariates for M -> Y
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
# Z[,2] <- (Z[,2]-mean(Z[,2]))/sd(Z[,2]) # scaled age, so we add a note that age is scaled before modeling
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
M_raw <- matrix(0,n,p) # the high-dimensional compositional mediator matrix 
for (i in 1:n){
  A1 <- M0[i,]*aX[i,]/sum(M0[i,]*aX[i,]) # the definition of  perturbation operator for M
  M_raw[i,] <- A1*U1[i,]/sum(A1*U1[i,]) # the definition of  perturbation operator for M
}
U2 <- rnorm(n,0,2) # the error term for Y in Eq. (17)
Y <- 1+ 0.5*trt + log(M_raw)%*%t(b) + Z%*%(eta)+ U2 # the outcome Y
OTU <- M_raw # high-dimensional compositional matrix

colnames(OTU) <- paste0("M", 1:50)
rownames(OTU) <- paste0("S", 1:200)

pheno <- data.frame(Treatment = trt, Outcome = Y, Sex = Z[,1], Age = Z[,2])

Example4 <- list(PhenoData = pheno, Mediator = OTU)

usethis::use_data(Example4, overwrite = TRUE)
