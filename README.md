#PGS: Penalized GEE with Grid Search

Penalized GEE model for association study of high-dimensional genomic data. 

##Download PGS R package (current version 0.2.0):

[PGS for Mac OS X](https://github.com/YinanZheng/PGS/releases/download/PGS_0.2.0/PGS_0.2.0.tgz)

[PGS for Windows (32/64)](https://github.com/YinanZheng/PGS/releases/download/PGS_0.2.0/PGS_0.2.0.zip)

[PGS for Linux (86/64)](https://github.com/YinanZheng/PGS/releases/download/PGS_0.2.0/PGS_0.2.0.tar.gz)

## Installation in R session

_`# First check and install dependencies:`_

    list.of.dependencies <- c("doParallel", "lme4", "geepack", "ggplot2", "reshape2")
    new.packages <- list.of.dependencies[!(list.of.dependencies %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages) else cat("Dependencies are ready!\n")
    
_`# Install PGS:`_

    install.packages("PGS_0.2.0.tgz", repo = NULL) # Mac OS X
    install.packages("PGS_0.2.0.zip", repo = NULL) # Windows
    install.packages("PGS_0.2.0.tar.gz", repo = NULL) # Linux

##Wiki & Examples:

[Wiki: PGS](https://github.com/YinanZheng/PGS/wiki)

[Example: Micro RNA expression and lung function](https://github.com/YinanZheng/PGS/wiki/Example:-miRNA-expression-and-lung-function)



##Citation:
1.	Zheng Y, Fei Z, Zhang W, Starren JB, Liu L, Baccarelli AA, Li Y, Hou L. PGS: a tool for association study of high-dimensional microRNA expression data with repeated measures. Bioinformatics. 2014;30(19):2802-7. doi: 10.1093/bioinformatics/btu396. PubMed PMID: 24947752; PubMed Central PMCID: PMC4173025. http://www.ncbi.nlm.nih.gov/pubmed/24947752

2.	Wang L, Zhou J, Qu A. Penalized generalized estimating equations for high-dimensional longitudinal data analysis. Biometrics. 2012;68(2):353-60. doi: 10.1111/j.1541-0420.2011.01678.x. PubMed PMID: 21955051. http://www.ncbi.nlm.nih.gov/pubmed/21955051

##Contact package maintainer:
Yinan Zheng 

Email: y-zheng@northwestern.edu



