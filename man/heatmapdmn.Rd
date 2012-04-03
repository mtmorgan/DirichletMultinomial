\name{heatmapdmn}
\alias{heatmapdmn}
\title{Heatmap representation of samples assigned to Dirichlet components.}
\description{
  Produce a heat map summarizing count data, grouped by Dirichlet component.
}
\usage{
heatmapdmn(count, fit1, fitN, ntaxa = 30, ...,
    transform = sqrt, lblwidth = 0.2 * nrow(count), col = .gradient)
}

\arguments{

  \item{count}{A matrix of sample x taxon counts, as supplied to
    \code{\link{dmn}}.}

  \item{fit1}{An instance of class \code{dmn}, from a model fit to a
    single Dirichlet component, \code{k=1} in \code{\link{dmn}}.}

  \item{fitN}{An instance of class \code{dmn}, from a model fit to
    \code{N != 1} components, \code{k=N} in \code{\link{dmn}}.}

  \item{ntaxa}{The \code{ntaxa} most numerous taxa to display counts for.}

  \item{\dots}{Additional arguments, ignored.}

  \item{transform}{Transformation to apply to count data prior to
    visualization; this does \emph{not} influence mixture membership or
    taxnomic ordering.}

  \item{lblwidth}{The proportion of the plot to dedicate to taxanomic
    labels, as a fraction of the number of samples to be plotted.}

  \item{col}{The colors used to display (possibly transformed, by
    \code{transform}) count data, as used by \code{\link{image}}.}

}

\details{

  Columns of the heat map correspond to samples. Samples are grouped by
  Dirichlet component, with average (Dirichlet) components summarized as
  a separate wide column. Rows correspond to taxonomic groups, ordered
  based on contribution to Dirichlet components.

}

\author{Martin Morgan \url{mailto:mtmorgan@fhcrc.org}}

\examples{
## counts
fl <- system.file(package="DirichletMultinomial", "extdata",
                  "Twins.csv")
count <- t(as.matrix(read.csv(fl, row.names=1)))

## all and best-fit clustering
data(fit)
lplc <- sapply(fit, laplace)
best <- fit[[which.min(lplc)]]

heatmapdmn(count, fit[[1]], best, 30)
}

\keyword{manip}
