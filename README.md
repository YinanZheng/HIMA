#HIMA: High-dimensional Mediation Analysis

HIMA is used to estimate and test high-dimensional mediation effects.

##Download HIMA R package (current version 0.0.1):

[HIMA for Mac OS X] coming soon

[HIMA for Windows (32/64)](https://github.com/YinanZheng/HIMA/releases/download/HIMA_v0.0.1/HIMA_0.0.1.zip)

[HIMA for Linux (86/64)] coming soon

## Installation in R session

_`# First check and install dependencies:`_

    list.of.dependencies <- c("glmnet","doParallel")
    new.packages <- list.of.dependencies[!(list.of.dependencies %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages) else cat("Dependencies are ready!\n")
    
_`# Install HIMA:`_

    install.packages("HIMA_0.0.1.tgz", repo = NULL) # Mac OS X
    install.packages("HIMA_0.0.1.zip", repo = NULL) # Windows
    install.packages("HIMA_0.0.1.tar.gz", repo = NULL) # Linux

##Wiki & Examples:

[Wiki: HIMA]

[Example: Simulation data]



##Citation:

Haixiang Zhang, Yinan Zheng, Zhou Zhang, Tao Gao, Grace Yoon, Wei Zhang, Lifang Hou, Andrea Baccarelli and Lei Liu. Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies. Bioinformatics 2016 (under review)

##Contact package maintainer:
Yinan Zheng 

Email: y-zheng@northwestern.edu



