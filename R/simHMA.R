#' Simulation Data Generator for High-dimensional Mediation Analyais
#' 
#' \code{simHMA} is used to simulate data for High-dimensional Mediation Analyais
#' 
#' @param y.vect a vector of dependent variable.
#' @param id.vect a vector of subjuect ID.
#' 
#' @seealso see \code{\link{hma}} to obtain proper ranked variables; see \code{\link{pgsfit.obj}} for class methods.
#' 
#' @examples
#' ### Dataset preview
#' BJdata()

n <- 100  # sample size
p <- 5000 # the dimension of covariates

a <- matrix(0,1,p) # the regression coefficients \alpha
b <- matrix(0,1,p) # the regression coefficients \beta
c <- runif(p,0,2)

a[1]  <- 0.5
a[2]  <- 0.5
a[3]  <- 0.3

##
b[1] <- 0.5
b[2] <- 1.2
b[3] <- 0.3
###  

M <- matrix(0,n,p) #  the mediators matrix with n by p
x <-  rnorm(n, mean=1, sd=1.5)
e <- matrix(0,1,n)
for (i in 1:p){
  e <- rnorm(n, 0, 1.2)
  y <- a[i]*x + c[i] + e
  M[,i] <- t(t(y))
}
##
X <- t(t(x)) #  the indepent covariates X with n by 1
X_M <- cbind(M,X) # the matrix combine M with X_2, 
##
E <- rnorm(n,0,1)
B <- c(b, 0.5) # (p+1) times 1
Y <-  1 + X_M%*%t(t(B)) + t(t(E)) #  the response with  n by 1


#  the following begin to compute the solution

a_est <- matrix(0,1,p) # the estimation for a

for (i in 1:p){
  random_data <- data.frame(y = M[,i], x1= X)
  solution <- lm(formula = y~x1, data = random_data)
  Est_coef <- coefficients(solution)
  est_hat <- as.numeric(Est_coef)
  a_est[i] <- est_hat[2]
}


##  the following is the ISI idea

d <- round(n/log10(n))  # the number of interested parameter

A_abs <- abs(a_est)

A_sort <- sort(A_abs)

E_ISI <- A_sort[p-d+1]

ID <- which(A_abs >= E_ISI) # the index of interested parameter


####
a_est_ISI <- a_est[ID]

M_ISI <- M[,ID]

X_ISI <- cbind(M_ISI,X)

fit = glmnet(X_ISI, Y, alpha = 0.2)
cv.fit <- cv.glmnet(X_ISI, Y, alpha = 0.2)
Coefficients <- coef(fit, s = cv.fit$lambda.min)
b_est_ISI <- Coefficients[2:(d+1)]  # the estimator for b



#### the following begin the Permutation

PN <- 1000 # the number of permutation
a_est_per <- matrix(0,PN,d)
b_est_per <- matrix(0,PN,d)
for (k in 1: PN){
  
  for (j in 1:d){
    random_data <- data.frame(y = sample(M_ISI[,j]), x1= X)
    solution <- lm(formula = y~x1, data = random_data)
    Est_coef <- coefficients(solution)
    est_hat <- as.numeric(Est_coef)
    a_est_per[k,j] <- est_hat[2]  # the estimaaiton for a_1 after permutation
  }
  ###
  
  fit = glmnet(X_ISI, sample(Y), alpha = 0.2)
  cv.fit <- cv.glmnet(X_ISI, sample(Y), alpha = 0.2)
  Coefficients <- coef(fit, s = cv.fit$lambda.min)
  b_est_per[k,] <- Coefficients[2:(d+1)]  # the estimator for b
} # the end of permutation




P_matrix   <- matrix(0,PN,d)
P_matrix_b <- matrix(0,PN,d)
for (i in 1:PN){
  P_matrix[i,]    <- as.numeric(abs(a_est_per[i,]) >= abs(a_est_ISI))
  P_matrix_b[i,]  <- as.numeric(abs(b_est_per[i,]) >= abs(b_est_ISI))
}


P_value_a <- colMeans(P_matrix)
P_value_b <- colMeans(P_matrix_b)

p_mix  <-  rbind(P_value_a,P_value_b)

P_value <- apply(p_mix,2,max) # the p-value

###### out-put results
print(ID)
##############
print(a_est_ISI)
##############
print(b_est_ISI)
##############
print(P_value)

