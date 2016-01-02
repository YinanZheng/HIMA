#Internal function: parallel computing check

checkParallel <- function(program.name, parallel, ncore, verbose = TRUE)
{
  if (parallel == TRUE & ncore > 1)
  { 
    if(ncore > detectCores())
    {
      cat(paste0("You requested ", ncore, " cores. There are only ", detectCores()," in your machine!"),'\n')
      ncore = detectCores()
    }
    if(verbose) cat(paste0("Start running ",program.name," with ", ncore, " cores in parallel...   (",Sys.time(),")\n"))
    if(getDoParWorkers() != ncore) registerDoParallel(ncore)
  } else {
    if(vebose) cat(paste0("Start running ",program.name," with single core...   (",Sys.time(),")\n"))
    registerDoSEQ()
  }
}
