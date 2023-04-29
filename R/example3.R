#' Example dataset 3 for HIMA: Survival data outcome
#'
#' A \code{list} dataset containing a phenotype dataset and a mediator
#' dataset (n=200 participants, p=200 biomarkers). The variables in the phenotype are as follows:
#'
#' \itemize{
#'   \item Treatment: treated (value = 1) or not treated (value = 0)
#'   \item Status: Status indicator: dead (value = 1) or alive (value = 0)
#'   \item Time: time to event
#'   \item Sex: female (value = 1) or male (value = 0)
#'   \item Age: Age of the participant
#' }
#'
#' The datasets are simulated using parameters generated from a real dataset. The code
#' used to generate the data can be found in /inst/script folder of the package.
#' 
#' @return A list containing two objects: \code{PhenoData} and \code{Mediator}
"Example3"
