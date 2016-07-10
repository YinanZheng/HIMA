## Internal function: doOne code generater

doOneGen <- function(model.text, colind.text){
  L = length(eval(parse(text = colind.text)))
  script = paste0("doOne <- function(i, datarun, Ydat){datarun$Mone <- Ydat[,i]; model <- ", model.text,";if('try-error' %in% class(model)) b <- rep(NA, ",L,") else { res=summary(model)$coefficients; b <- res[2,",colind.text,"]};invisible(b)}")
  return(script)
}
