## HIMA: High-dimensional Mediation Analysis

![GitHub release](https://img.shields.io/badge/release-v2.2.0-blue.svg)

*HIMA* is an R package for estimating and testing high-dimensional mediation effects in omic studies.

## Installation 

HIMA is now available in R CRAN repository
```r
## Install HIMA
install.packages("HIMA")
```

To install from GitHub
```r
## Install HIMA from GitHub
library(devtools)
install_github("yinanzheng/HIMA")
```

If package "qvalue" is not found, please first install "qvalue" package through Bioconductor: https://www.bioconductor.org/packages/release/bioc/html/qvalue.html
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")
```

## Citation:

Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J, Just A, Colicino E, Vokonas P, Zhao L, Lv J, Baccarelli A, Hou L, Liu L. Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies. Bioinformatics. 2016;32(20):3150-4. DOI: 10.1093/bioinformatics/btw351. PubMed PMID: 27357171; PMCID: PMC5048064.

Zhang H, Zheng Y, Hou L, Liu L. Mediation Analysis for Survival Data with High-Dimensional Mediators. Bioinformatics. 2021;37(21):3815-21. DOI: 10.1093/bioinformatics/btab564. PubMed PMID: 34343267; PMCID: PMC8570823.

Zhang H, Chen J, Feng Y, Wang C, Li H, Liu L. Mediation effect selection in high-dimensional and compositional microbiome data. Stat Med. 2021;40(4):885-96. DOI: 10.1002/sim.8808. PubMed PMID: 33205470; PMCID: PMC7855955.

## Contact package maintainer:

Yinan Zheng 

Email: y-zheng@northwestern.edu
