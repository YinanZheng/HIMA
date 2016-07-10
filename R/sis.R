#Internal function: Sure Independent Screening

himasis <- function(Y, M, X, COV, modelstatement, scale = TRUE, parallel=TRUE, ncore = 2)
{
  L.M = ncol(M); M.names = colnames(M)
  
  if(scale)
    M_scale = scaleto(M)
  else M_scale = M
  
  X <- data.frame(X)
  X <- data.frame(model.matrix(~.,X))[,-1]
  

  if(is.null(COV)){
    datarun = data.frame(Y = Y, Mone = NA, X = X)
    modelstatement = modelstatement
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~.,COV))[,-1]
    conf.names = colnames(COV)   
    datarun = data.frame(Y = Y, Mone = NA, X = X, COV = COV)
    modelstatement =  eval(parse(text=(paste0(modelstatement, "+", paste0(paste0("COV.",conf.names),collapse = "+")))))
  }

  doOne <- eval(parse(text = doOneGen("try(lm(modelstatement, data = datarun))","c(1,4)")))
    
  checkParallel("Sure Independent Screening", parallel, ncore, verbose = F)
  
  results <- foreach(n = idiv(L.M, chunks = ncore), M_chunk = iblkcol_lag(M, chunks = ncore),.combine = 'cbind') %dopar% {
    sapply( seq_len(n), doOne, datarun, M_chunk )
  }
  
  colnames(results) = M.names
  return(results)
}