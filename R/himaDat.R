#' Continuous Outcome Dataset for HIMA Demo
#'
#' A dataset containing phenotype data and high-dimensional mediators for continuous outcome analysis.
#' The dataset was simulated using parameters generated from real data.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{PhenoData}{A data frame containing:
#'     \describe{
#'       \item{Treatment}{treated (value = 1) or not treated (value = 0).}
#'       \item{Outcome}{a normally distributed continuous outcome variable.}
#'       \item{Sex}{female (value = 1) or male (value = 0).}
#'       \item{Age}{age of the participant.}
#'     }
#'   }
#'   \item{Mediator}{A matrix of high-dimensional mediators (rows: samples, columns: variables).}
#' }
#'
#' @examples
#' data(ContinuousOutcome)
#' head(ContinuousOutcome$PhenoData)
#'
"ContinuousOutcome"

#' Binary Outcome Dataset for HIMA Demo
#'
#' A dataset containing phenotype data and high-dimensional mediators for binary outcome analysis.
#' The dataset was simulated using parameters generated from real data.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{PhenoData}{A data frame containing:
#'     \describe{
#'       \item{Treatment}{treated (value = 1) or not treated (value = 0).}
#'       \item{Disease}{binary outcome: diseased (value = 1) or healthy (value = 0).}
#'       \item{Sex}{female (value = 1) or male (value = 0).}
#'       \item{Age}{age of the participant.}
#'     }
#'   }
#'   \item{Mediator}{A matrix of high-dimensional mediators (rows: samples, columns: variables).}
#' }
#'
#' @examples
#' data(BinaryOutcome)
#' head(BinaryOutcome$PhenoData)
#'
"BinaryOutcome"

#' Survival Outcome Dataset for HIMA Demo
#'
#' A dataset containing phenotype data and high-dimensional mediators for survival outcome analysis.
#' The dataset was simulated using parameters generated from real data.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{PhenoData}{A data frame containing:
#'     \describe{
#'       \item{Treatment}{treated (value = 1) or not treated (value = 0).}
#'       \item{Status}{status indicator: dead (value = 1) or alive (value = 0).}
#'       \item{Time}{time to the event or censoring.}
#'       \item{Sex}{female (value = 1) or male (value = 0).}
#'       \item{Age}{age of the participant.}
#'     }
#'   }
#'   \item{Mediator}{A matrix of high-dimensional mediators (rows: samples, columns: variables).}
#' }
#'
#' @examples
#' data(SurvivalData)
#' head(SurvivalData$PhenoData)
#'
"SurvivalData"

#' Compositional Mediator Dataset for HIMA Demo
#'
#' A dataset containing phenotype data and high-dimensional compositional mediators (e.g., microbiome).
#' The dataset was simulated using parameters generated from real data.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{PhenoData}{A data frame containing:
#'     \describe{
#'       \item{Treatment}{treated (value = 1) or not treated (value = 0).}
#'       \item{Outcome}{a normally distributed continuous outcome variable.}
#'       \item{Sex}{female (value = 1) or male (value = 0).}
#'       \item{Age}{age of the participant.}
#'     }
#'   }
#'   \item{Mediator}{A matrix of high-dimensional compositional mediators (rows: samples, columns: variables).}
#' }
#'
#' @examples
#' data(MicrobiomeData)
#' head(MicrobiomeData$PhenoData)
#'
"MicrobiomeData"

#' Quantile Mediation Dataset for HIMA Demo
#'
#' A dataset containing phenotype data and high-dimensional mediators for quantile mediation analysis.
#' The dataset was simulated using parameters generated from real data.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{PhenoData}{A data frame containing:
#'     \describe{
#'       \item{Treatment}{treated (value = 1) or not treated (value = 0).}
#'       \item{Outcome}{an abnormally distributed continuous outcome variable.}
#'       \item{Sex}{female (value = 1) or male (value = 0).}
#'       \item{Age}{age of the participant.}
#'     }
#'   }
#'   \item{Mediator}{A matrix of high-dimensional mediators (rows: samples, columns: variables).}
#' }
#'
#' @examples
#' data(QuantileData)
#' head(QuantileData$PhenoData)
#'
"QuantileData"
