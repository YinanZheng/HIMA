#' Example dataset 1 for HIMA: Continuous outcome 
#'
#' A \code{list} dataset containing a phenotype dataset and a high-dimension mediator
#' dataset (n=300 participants, p=300 biomarkers). The variables in the phenotype are as follows:
#'
#' \itemize{
#'   \item Treatment: treated (value = 1) or not treated (value = 0)
#'   \item Outcome: outcome of the treatment- a normally distributed continuous variable
#'   \item Sex: female (value = 1) or male (value = 0)
#'   \item Age: Age of the participant
#' }
#'
#' The datasets are simulated using parameters generated from a real dataset. The code
#' used to generate the data can be found in /inst/script folder of the package.
#'
#' @return A list containing two objects: \code{PhenoData} and \code{Mediator}
"Example1"
