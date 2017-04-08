## HIMA: High-dimensional Mediation Analysis
[![GitHub release](https://img.shields.io/badge/release-v1.0.4-blue.svg)](https://github.com/YinanZheng/HIMA/releases)

*HIMA* is an [R](http://en.wikipedia.org/wiki/R_%28programming_language%29) package for estimating and testing high-dimensional mediation effects in genomic/epigenomic studies.

## Installation 

In R console,
```r
## Install REMP
library(devtools)
install_github("YinanZheng/HIMA",
               dependencies=TRUE)
               
## If SSL cert verification failure
library(RCurl)
library(httr)
set_config( config( ssl_verifypeer = 0L ) )
install_github("YinanZheng/HIMA",
               dependencies=TRUE)
```

##Citation:

Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J, Just A, Colicino E, Vokonas P, Zhao L, Lv J, Baccarelli A, Hou L, Liu L. Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies. Bioinformatics. 2016. DOI: [http://dx.doi.org/10.1093/bioinformatics/btw351](http://dx.doi.org/10.1093/bioinformatics/btw351). [PubMed PMID: 27357171](http://www.ncbi.nlm.nih.gov/pubmed/?term=27357171).

##Contact package maintainer:
Yinan Zheng 

Email: y-zheng@northwestern.edu
