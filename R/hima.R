hima <- function(...)
{
  
  a_est <- matrix(0,1,p) # the estimation for a
  
  M = simdat$M
  X = simdat$X
  Y = simdat$Y
  n = simdat$n
  p = simdat$p
  
  sis.res = sis(X, id.vect = 1:n, M)
  id = as.character(coef(sis.res)[1:d,]$name)
  sort(as.numeric(substring(id, 2, nchar(id))))
  
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
  
  A_sort <- sort(A_abs, decreasing = T)
  
  E_ISI <- A_sort[d]
  
  ID <- which(A_abs >= E_ISI) # the index of interested parameter
  
  head(ID)
  
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
  
  
  
  
  P_matrix_a   <- matrix(0,PN,d)
  P_matrix_b <- matrix(0,PN,d)
  for (i in 1:PN){
    P_matrix_a[i,]    <- as.numeric(abs(a_est_per[i,]) >= abs(a_est_ISI))
    P_matrix_b[i,]  <- as.numeric(abs(b_est_per[i,]) >= abs(b_est_ISI))
  }
  
  
  P_value_a <- colMeans(P_matrix_a)
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
  
  
  
}