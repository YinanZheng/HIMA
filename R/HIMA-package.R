#' @title High-dimensional Mediation Analysis
#'
#' @description
#' HIMA is an R package for estimating and testing high-dimensional mediation effects in omic studies.
#' 
#' \tabular{ll}{ Package: \tab HIMA\cr Type: \tab Package\cr
#' Version: \tab 1.1.0\cr Date: \tab 2021-05-12\cr License: \tab GPL-3\cr }
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
#' @references 1. Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J, Just A, Colicino E, Vokonas P, Zhao L, 
#' Lv J, Baccarelli A, Hou L, Liu L. Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies. 
#' Bioinformatics. 2016. DOI: \url{http://dx.doi.org/10.1093/bioinformatics/btw351}. PubMed PMID: 27357171.
#' 
#' 2. Zhang H, Zheng Y, Hou L, Liu L. Mediation Analysis for Survival Data with High-Dimensional Mediators. 
#' Bioinformatics. 2021 (under review).
#' 
#' @keywords package
#' 
#' @import ncvreg
#' @import glmnet
#' 
#' @importFrom utils packageVersion
#' @importFrom stats p.adjust rnorm runif rbinom glm model.matrix BIC coef ls.diag lsfit pnorm
#' @importFrom iterators idiv nextElem
#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar% registerDoSEQ getDoParWorkers
#' @importFrom doParallel registerDoParallel
#' @importFrom HDMT null_estimation fdr_est
#' @importFrom survival Surv coxph
#' @importFrom MASS mvrnorm
#' 
NULL