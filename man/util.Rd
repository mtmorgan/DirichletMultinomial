\name{Utilities}
\alias{csubset}
\title{Helpful utility functions}

\description{
  \code{csubset} creates a subset of a count matrix, based on identity
  of column phenotypes to a specified value.
}
\usage{
csubset(val, x, pheno, cidx = TRUE)
}

\arguments{

  \item{val}{\code{character(1)} specifying the subset of
    \code{phenotype} to select.}

  \item{x}{A matrix of counts, with rows corresponding to samples and
    columns to taxonomic groups.}

  \item{pheno}{A \code{character()} vector of length equal to the number
    of rows in \code{count}, indicating the phenotype of the
    corresponding sample.}

  \item{cidx}{A \code{logical(1)} indicating whether columns (taxa) with
    zero counts in the count matrix following removal of taxa not
    satisfying \code{pheno \%in\% val} should be removed. 
    \code{cidx=FALSE} removes the 0-count columns.}

}

\value{

  A \code{matrix} of counts, with rows satisfying \code{pheno \%in\%
  val} and with columns equal either to \code{ncol(x)} (when
  \code{cidx=TRUE}) or the number of columns with non-zero counts after
  row subsetting (\code{cidx=FALSE}).

}
\author{Martin Morgan \url{mailto:mtmorgan@fhcrc.org}}

\examples{
## count matrix
fl <- system.file(package="DirichletMultinomial", "extdata",
                  "Twins.csv")
count <- t(as.matrix(read.csv(fl, row.names=1)))

## phenotype
fl <- system.file(package="DirichletMultinomial", "extdata",
                  "TwinStudy.t")
pheno0 <- scan("~/DirichletMultinomial/inst/extdata/TwinStudy.t")
lvls <- c("Lean", "Obese", "Overwt")
pheno <- factor(lvls[pheno0 + 1], levels=lvls)
names(pheno) <- rownames(count)

## subset
dim(count)
sum("Lean" == pheno)
dim(csubset("Lean", count, pheno))
dim(csubset("Lean", count, pheno, cidx=FALSE))
}

\keyword{manip}
