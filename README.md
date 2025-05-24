## HIMA: High-dimensional Mediation Analysis

![GitHub release](https://img.shields.io/badge/release-v2.3.2-blue.svg)

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
## Documentation

Explore the tutorial for the `HIMA` package, which provides detailed usage examples and guidelines, at the following tutorial link:

[HIMA Tutorial](https://YinanZheng.github.io/HIMA/articles/hima-vignette.html)

## Citation:

1. Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J, Just A, Colicino E, Vokonas P, Zhao L, 
Lv J, Baccarelli A, Hou L, Liu L. Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies. 
Bioinformatics. 2016. DOI: 10.1093/bioinformatics/btw351. PMID: 27357171; PMCID: PMC5048064

2. Zhang H, Zheng Y, Hou L, Zheng C, Liu L. Mediation Analysis for Survival Data with High-Dimensional Mediators. 
Bioinformatics. 2021. DOI: 10.1093/bioinformatics/btab564. PMID: 34343267; PMCID: PMC8570823

3. Zhang H, Chen J, Feng Y, Wang C, Li H, Liu L. Mediation effect selection in high-dimensional and compositional microbiome data. 
Stat Med. 2021. DOI: 10.1002/sim.8808. PMID: 33205470; PMCID: PMC7855955

4. Zhang H, Chen J, Li Z, Liu L. Testing for mediation effect with application to human microbiome data. 
Stat Biosci. 2021. DOI: 10.1007/s12561-019-09253-3. PMID: 34093887; PMCID: PMC8177450

5. Perera C, Zhang H, Zheng Y, Hou L, Qu A, Zheng C, Xie K, Liu L. HIMA2: high-dimensional mediation analysis and its application in epigenome-wide DNA methylation data. 
BMC Bioinformatics. 2022. DOI: 10.1186/s12859-022-04748-1. PMID: 35879655; PMCID: PMC9310002

6. Zhang H, Hong X, Zheng Y, Hou L, Zheng C, Wang X, Liu L. High-Dimensional Quantile Mediation Analysis with Application to a Birth 
Cohort Study of Mother–Newborn Pairs. Bioinformatics. 2024. DOI: 10.1093/bioinformatics/btae055. PMID: 38290773; PMCID: PMC10873903

7. Bai X, Zheng Y, Hou L, Zheng C, Liu L, Zhang H. An Efficient Testing Procedure for High-dimensional Mediators with FDR Control. 
Statistics in Biosciences. 2024. DOI: 10.1007/s12561-024-09447-4.

## Contact package maintainer:

Yinan Zheng 

Email: y-zheng@northwestern.edu
