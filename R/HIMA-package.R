#' @title High-dimensional Mediation Analysis
#'
#' @description
#' HIMA is an R package for estimating and testing high-dimensional mediation effects in genomic/epigenomic studies.
#' 
#' \tabular{ll}{ Package: \tab HIMA\cr Type: \tab Package\cr
#' Version: \tab 1.0.6\cr Date: \tab 2018-02-28\cr License: \tab GPL-3\cr }
#' 
#' @name HIMA-package
#' 
#' @aliases HIMA-package HIMA
#' 
#' @docType package
#' 
#' @author 
#' Yinan Zheng \email{y-zheng@@northwestern.edu},
#' Haixiang Zhang \email{haixiang.zhang@@tju.edu.cn},  
#' Lifang Hou \email{l-hou@@northwestern.edu}
#' Lei liu \email{lei.liu@@wustl.edu}
#' 
#' Maintainer: Yinan Zheng \email{y-zheng@@northwestern.edu}
#' 
#' @references Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J, Just A, Colicino E, Vokonas P, Zhao L, 
#' Lv J, Baccarelli A, Hou L, Liu L. Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies. 
#' Bioinformatics. 2016. DOI: \url{http://dx.doi.org/10.1093/bioinformatics/btw351}. PubMed PMID: 27357171.
#' 
#' @keywords package
#' 
#' @import ncvreg
#' 
#' @importFrom stats p.adjust rnorm runif rbinom glm model.matrix BIC coef
#' @importFrom iterators idiv nextElem
#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar% registerDoSEQ getDoParWorkers
#' @importFrom doParallel registerDoParallel

NULL