#HMA: High-dimensional Mediation Analysis

High-dimensional Mediation Analysis

##Download HMA R package (current version 0.0.1):

[HMA for Mac OS X](https://github.com/YinanZheng/HMA/releases/download/HMA_0.0.1/HMA_0.0.1.tgz)

[HMA for Windows (32/64)](https://github.com/YinanZheng/HMA/releases/download/HMA_0.0.1/HMA_0.0.1.zip)

[HMA for Linux (86/64)](https://github.com/YinanZheng/HMA/releases/download/HMA_0.0.1/HMA_0.0.1.tar.gz)

## Installation in R session

_`# First check and install dependencies:`_

    list.of.dependencies <- c("glmnet")
    new.packages <- list.of.dependencies[!(list.of.dependencies %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages) else cat("Dependencies are ready!\n")
    
_`# Install HMA:`_

    install.packages("HMA_0.0.1.tgz", repo = NULL) # Mac OS X
    install.packages("HMA_0.0.1.zip", repo = NULL) # Windows
    install.packages("HMA_0.0.1.tar.gz", repo = NULL) # Linux

##Wiki & Examples:

[Wiki: HMA]

[Example: Simulation data]



##Citation:


##Contact package maintainer:
Yinan Zheng 

Email: y-zheng@northwestern.edu



