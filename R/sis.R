#Internal function: Sure Independent Screening

sis <- function(Y, x.vect, COV=NULL, parallel=TRUE, ncore = 2)
{
  L.Y = ncol(Y); Y.names = colnames(Y)
  res.COV = NULL
  
  if (is.null(COV))
  { 
    datarun = data.frame(Yone = NA, x.vect=x.vect)
    modelstatement_MLR = modelstatement_GEE = Yone ~ x.vect
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~.,COV))[,-1]
    conf.names = colnames(COV)          
    datarun = data.frame(Yone = NA, x.vect=x.vect, COV = COV)
    modelstatement_MLR = eval(parse(text=(paste0("Yone ~ x.vect +", paste0(paste0("COV.",conf.names),collapse = "+")))))
  }
  
  modelstatement = modelstatement_MLR
  doOne <- eval(parse(text = doOneGen("try(lm(modelstatement, data = datarun))","1")))
    
  checkParallel("Sure Independent Screening", parallel, ncore, verbose = F)
  
  results <- foreach(n = idiv(L.Y, chunks = ncore), Y_chunk = iblkcol_lag(Y, chunks = ncore),.combine = 'c') %dopar% {
    do.call('rbind',lapply( seq_len(n), doOne, datarun, Y_chunk) )
  }
  
  names(results) = Y.names
  return(results)
}