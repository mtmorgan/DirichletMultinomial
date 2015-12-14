\name{cvdmngroup}
\alias{cvdmngroup}

\title{
  Cross-validation on Dirichlet-Multinomial classifiers.
}
\description{

  Run cross-validation on Dirichlet-Multinomial generative classifiers.

}
\usage{
cvdmngroup(ncv, count, k, z, ..., verbose = FALSE,
    .lapply = parallel::mclapply)
}

\arguments{
  \item{ncv}{\code{integer(1)} number of cross-validation groups,
    between 2 and \code{nrow(count)}.}

  \item{count}{\code{matrix} of sample x taxon counts, subsets of which
    are used for training and cross-validation.}

  \item{k}{named \code{integer()} vector of groups and number of
    Dirichlet components; e.g., \code{c(Lean=1, Obese=3)} performs
    cross-validation for models with \code{k=1} Dirichlet components for
    the \sQuote{Lean} group, \code{k=3} Dirichlet components for
    \sQuote{Obese}.}

  \item{z}{True group assignment.}

  \item{\dots}{Additional arguments, passed to \code{\link{dmn}} during
    each cross-validation.}

  \item{verbose}{\code{logical(1)} indicating whether progress should be
    reported}

  \item{.lapply}{A function used to perform the outer cross-vaildation
    loop, e.g., \code{lapply} for calculation on a single processor,
    \code{parallel::mclapply} for parallel evaluation.}

}
\value{

  A \code{data.frame} summarizing classifications of test samples in
  cross-validation groups. Columns are:
  
  \item{group}{The cross-validation group in which the indivdual was
    used for testing.}

  \item{additional columns}{Named after classification groups, giving
    the posterior probability of assignment.}

}
\author{Martin Morgan \url{mailto:mtmorgan@fhcrc.org}}

\seealso{
  \code{\link{dmn}}, \link{DirichletMultinomial-package},
  \code{vignette("DirichletMultinomial")}
}
\examples{
data(xval)  ## result of following commands
head(xval)

\dontrun{
## count matrix
fl <- system.file(package="DirichletMultinomial", "extdata",
                  "Twins.csv")
count <- t(as.matrix(read.csv(fl, row.names=1)))

## phenotype
fl <- system.file(package="DirichletMultinomial", "extdata",
                  "TwinStudy.t")
pheno0 <- scan(fl)
lvls <- c("Lean", "Obese", "Overwt")
pheno <- factor(lvls[pheno0 + 1], levels=lvls)
names(pheno) <- rownames(count)

## subset
keep <- c("Lean", "Obese")
count <- count[pheno %in% keep,]
pheno <- factor(pheno[pheno %in% keep], levels=keep)

## cross-validation, single Dirichlet component for Lean, 3 for Obese
xval <- cvdmngroup(nrow(count), count, c(Lean=1, Obese=3), pheno,
                   verbose=TRUE, mc.preschedule=FALSE)
}}

\keyword{stats}
