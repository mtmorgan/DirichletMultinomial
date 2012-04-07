\name{DMN-class}
\Rdversion{1.1}
\docType{class}
\alias{DMN-class}

\title{Class \code{"DMN"}}
\description{Result from fitting a Dirichlet-Multinomial model.}
\section{Objects from the Class}{
  Objects can be created by calls to \code{\link{dmn}}..
}
\section{Slots}{

  The contents of a slot is usually retrieved via the methods described
  on the \code{\link{mixture}} help page.

  \describe{

    \item{goodnessOfFit}{NLE, LogDet, Laplace, AIC, and BIC criteria
      assessing goodness-of-fit.}

    \item{group}{\code{matrix} of dimension samples x \code{k},
      providing the Dirichlet parameter vectors.}

    \item{mixture}{\describe{
	
	\item{Weight}{\code{numeric()} of length \code{k}, with relative
	  weight of each component.}
    }}

    \item{fit}{\describe{

	\item{Lower}{\code{matrix()} of dimension taxa x \code{k} with
	  95\% lower bounds on Dirichlet component vector estimates.}
	
	\item{Estimate}{\code{matrix()} of dimension taxa x \code{k} with
	  Dirichlet component vector estimates.}

	\item{Upper}{\code{matrix()} of dimension taxa x \code{k} with
	  95\% upper bounds on Dirichlet component vector estimates.}

    }}

  }
}
\section{Methods}{

  See the \code{\link{mixture}} help page.

}
\author{
  Martin Morgan \url{mailto:mtmorgan@fhcrc.org}
}

\seealso{
  \code{\link{dmn}}, \code{\link{mixture}}.
}
\examples{
data(fit)
fit[[4]]
}
\keyword{classes}
