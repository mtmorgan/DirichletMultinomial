\name{dmngroup}
\alias{dmngroup}

\title{
  Dirichlet-Multinomial generative classifiers.
}
\description{

  Fit Dirichlet-Multinomial generative classifiers to groups (rows)
  within a sample x taxon count matrix.

}
\usage{

dmngroup(count, group, k, ..., simplify = TRUE,
    .lapply = parallel::mclapply)

}
\arguments{

  \item{count}{\code{matrix()} of sample x taxon counts.}

  \item{group}{\code{factor()} or vector to be coerced to a factor, with
    as many elements as there are rows in \code{count}, indicating the
    group to which the corresponding sample belongs.}

  \item{k}{\code{integer()}, the number(s) of Dirichlet components to
    fit.}

  \item{\dots}{Additional arguments, passed to \code{\link{dmn}}.}

  \item{simplify}{Return only the best-fit model for each group?}

  \item{.lapply}{An \code{lapply}-like function for application of group
    x k fits.}

}
\details{

  This function divided \code{count} into groups defined by
  \code{group}, creates all combinations of \code{group} x \code{k},
  and evaluates each using \code{\link{dmn}}. When \code{simplify=TRUE},
  the best (Laplace) fit is selected for each group.

}
\value{
  An object of class \code{dmngroup}, a list of fitted models of class
  \code{\link{dmn}}. When \code{simplify=TRUE}, elements are named by
  the group to which they correspond.
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
  \code{\link{dmn}}, \link{DirichletMultinomial-package},
  \code{vignette("DirichletMultinomial")}
}
\examples{
## best fit for groups 'Lean' and 'Obese'; full example in vignette.
\dontrun{bestgrp <- dmngroup(count, pheno, k=1:5, verbose=TRUE, 
                    mc.preschedule=FALSE)
}
data(bestgrp)
bestgrp
bestgrp[["Obese"]]
}

\keyword{ manip }
