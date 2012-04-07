\name{data}
\alias{fit}
\alias{bestgrp}
\alias{xval}
\docType{data}
\title{Data objects used for examples and the vignette}
\description{

  These data objects correspond to steps in a typical work flow, as
  described in the vignette to this package. \code{fit} corresponds to
  \code{dmn} fits to different values of \code{k}. \code{bestgroup} is
  the result of the two-group generative classifier. \code{xval}
  summarizes leave-one-out cross validation of the classifier.

}
\usage{
data(fit)
data(bestgrp)
data(xval)
}
\format{

  \code{fit} is a list of seven \code{\linkS4class{DMN}} objects.

  \code{bestgrp} is a \code{\linkS4class{DMNGroup}} object.

  \code{xval} is a \code{data.frame} with columns corresponding to the
  cross-validation group membership and the Lean and Obese posterior
  probabilities.

}
\examples{
data(fit); fit[1:2]
plot(sapply(fit, laplace), type="b")
data(bestgrp); bestgrp
data(xval); head(xval, 3)
}
\keyword{datasets}
