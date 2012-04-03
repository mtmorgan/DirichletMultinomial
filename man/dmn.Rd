\name{dmn}
\alias{dmn}

\title{
  Fit Dirichlet-Multinomial models to count data.
}
\description{
  Fit Dirichlet-Multinomial models to a sample x taxon count matrix.
}
\usage{
dmn(count, k, verbose = FALSE, seed = runif(1, 0, .Machine$integer.max))
}
\arguments{
  \item{count}{\code{matrix()} of sample x taxon counts.}

  \item{k}{\code{integer(1)}, the number of Dirichlet components to fit.}

  \item{verbose}{\code{logical(1)} indicating whether progress in fit
    should be reported.}

  \item{seed}{\code{numeric(1)} random number seed.}

}
\details{
  This implements Dirichlet-multinomial mixture models describe in the
  package help page, \link{DirichletMultinomial-package}.
}
\value{
  An object of class \code{dmn}, with elements (elements are usually
  retrieved via functions defined in the package, not directly).

  \item{GoodnessOfFit}{NLE, LogDet, Laplace, AIC, and BIC criteria
    assessing goodness-of-fit.}

  \item{Group}{\code{matrix} of dimension samples x \code{k},
  providing the Dirichlet parameter vectors.}

  \item{Mixture}{\describe{
      
      \item{Weight}{\code{numeric()} of length \code{k}, with relative
	weight of each component.}
    }}

    \item{Fit}{\describe{

	\item{Lower}{\code{matrix()} of dimension taxa x \code{k} with
	  95\% lower bounds on Dirichlet component vector estimates.}
	
	\item{Estimate}{\code{matrix()} of dimension taxa x \code{k} with
	  Dirichlet component vector estimates.}

	\item{Upper}{\code{matrix()} of dimension taxa x \code{k} with
	  95\% upper bounds on Dirichlet component vector estimates.}

      }}
}

\references{

  Holmes I, Harris K, Quince C, 2012 Dirichlet Multinomial Mixtures:
  Generative Models for Microbial Metagenomics. PLoS ONE 7(2):
  e30126. doi:10.1371/journal.pone.0030126.

}
\author{
  Martin Morgan \url{mailto:mtmorgan@fhcrc.org}
}

\seealso{
  \link{DirichletMultinomial-package},
  \code{vignette("DirichletMultinomial")}
}
\examples{
data(fit)
## k = 1:7; full example in vignette
lplc <- sapply(fit, laplace)
plot(lplc, type="b")
fit[[which.min(lplc)]]
}
\keyword{ manip }
