% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/himaDat.R
\docType{data}
\name{himaDat}
\alias{himaDat}
\title{Example datasets for HIMA demo}
\format{
An object of class \code{list} of length 5.
}
\usage{
himaDat
}
\value{
A list of example datasets for HIMA demo and testing.\cr\cr
-->>
\bold{Example dataset 1 for HIMA: Continuous outcome}
<<--
\describe{
  \item{Treatment:}{treated (value = 1) or not treated (value = 0)}
  \item{Outcome:}{outcome of the treatment- a normally distributed continuous variable}
  \item{Sex:}{female (value = 1) or male (value = 0)}
  \item{Age:}{Age of the participant}
}
-->>
\bold{Example dataset 2 for HIMA: Binary outcome }
<<--
\describe{
  \item{Treatment:}{treated (value = 1) or not treated (value = 0)}
  \item{Disease:}{diseased (value = 1) or healthy (value = 0)}
  \item{Sex:}{female (value = 1) or male (value = 0)}
  \item{Age:}{Age of the participant}
}
-->>
\bold{Example dataset 3 for HIMA: Survival data outcome}
<<--
\describe{
  \item{Treatment:}{treated (value = 1) or not treated (value = 0)}
  \item{Status:}{Status indicator: dead (value = 1) or alive (value = 0)}
  \item{Time:}{time to event}
  \item{Sex:}{female (value = 1) or male (value = 0)}
  \item{Age:}{Age of the participant}
}
-->>
\bold{Example dataset 4 for HIMA: Compositional mediator (e.g., microbiome)}
<<--
\describe{
  \item{Treatment:}{treated (value = 1) or not treated (value = 0)}
  \item{Outcome:}{outcome of the treatment- a normally distributed continuous variable}
  \item{Sex:}{female (value = 1) or male (value = 0)}
  \item{Age:}{Age of the participant}
}
-->>
\bold{Example dataset 5 for HIMA: High-dimensional quantile mediation analysis}
<<--
\describe{
  \item{Treatment:}{treated (value = 1) or not treated (value = 0)}
  \item{Outcome:}{outcome of the treatment- abnormally distributed continuous variable}
  \item{Sex:}{female (value = 1) or male (value = 0)}
  \item{Age:}{Age of the participant}
}
}
\description{
A \code{list} containing datasets for various scenarios of HIMA. Each dataset contains a phenotype 
data frame and a high-dimension mediator data matrix. The datasets are simulated using parameters generated 
from real datasets. The code used to generate the data can be found in "/inst/script" folder of the source package.
}
\keyword{datasets}
