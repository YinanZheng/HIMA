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

#' Longitudinal Mediators with Survival Outcome Dataset for HIMA Demo
#'
#' A simulated dataset for demonstrating high-dimensional and longitudinal
#' mediation analysis with survival outcomes in a counting-process framework.
#' The data were generated under a longitudinal mediator model and a
#' piecewise-constant Weibull survival model, mimicking real-world analysis
#' settings.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{PhenoData}{A data frame where each row represents one interval
#'   (\code{tstart}, \code{tstop}) for a subject in counting-process format.
#'   It contains:
#'     \describe{
#'       \item{ID}{Subject identifier (may appear multiple times due to interval splitting).}
#'       \item{Tstart}{Start time of the interval.}
#'       \item{Tstop}{Stop time of the interval (event or censoring time).}
#'       \item{Status}{Event indicator for the interval (\code{1 = event}, \code{0 = no event}).}
#'       \item{Treatment}{Exposure variable for each subject.}
#'       \item{Sex}{Binary covariate: \code{1 = male}, \code{0 = female}.}
#'       \item{Age}{Age of the subject in years.}
#'     }
#'   }
#'   \item{Mediator}{A numeric matrix of high-dimensional longitudinal mediators
#'   aligned with the rows of \code{PhenoData}.  
#'   Columns correspond to mediator variables (\code{M1}, \code{M2}, …), and rows
#'   correspond to observation intervals in the counting-process setup.}
#' }
#'
#' @examples
#' data(SurvivalLongData)
#' head(SurvivalLongData$PhenoData)
#'
"SurvivalLongData"

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
