#' @title High-Dimensional Mediation Analysis for 'Omic' Data
#'
#' @description
#' HIMA is an R package for estimating and testing high-dimensional mediation effects in omic studies.
#' HIMA can perform high-dimensional mediation analysis on a wide range of omic data types as
#' potential mediators, including epigenetics, transcriptomics, proteomics, and metabolomics using function \code{hima} and 
#' microbiome data (function \code{microHIMA}). HIMA can also handle survival data (function \code{survHIMA}).
#' 
#' \tabular{ll}{ Package: \tab HIMA\cr Type: \tab Package\cr
#' Version: \tab 2.0.0\cr Date: \tab 2022-02-02\cr License: \tab GPL-3\cr }
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
#' Bioinformatics. 2016. DOI: 10.1093/bioinformatics/btw351. PMID: 27357171. PMCID: PMC5048064
#' 
#' 2. Zhang H, Zheng Y, Hou L, Zheng C, Liu L. Mediation Analysis for Survival Data with High-Dimensional Mediators. 
#' Bioinformatics. 2021. DOI: 10.1093/bioinformatics/btab564. PMID: 34343267. PMCID: PMC8570823
#' 
#' @keywords package
#' 
#' @import ncvreg
#' @import glmnet
#' 
#' @importFrom utils packageVersion
#' @importFrom stats p.adjust rnorm runif rbinom lm glm model.matrix BIC coef ls.diag lsfit pnorm complete.cases rgamma
#' @importFrom iterators idiv nextElem
#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar% registerDoSEQ getDoParWorkers
#' @importFrom doParallel registerDoParallel
#' @importFrom HDMT fdr_est
#' @importFrom survival Surv coxph
#' @importFrom MASS mvrnorm glm.nb
#' @importFrom hommel hommel discoveries
#' 
NULL