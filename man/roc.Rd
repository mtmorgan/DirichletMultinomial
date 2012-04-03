\name{roc}
\alias{roc}
\title{Summarize receiver-operator characteristics}
\description{

  Returns a \code{data.frame} summarizing the cummulative true- and
  false-positive probabilities from expected and observed
  classifications.

}

\usage{
roc(exp, obs, ...)
}

\arguments{

  \item{exp}{\code{logical()} vector of expected classifications to a
  particular group.}

  \item{obs}{Predicted probability of assignment to the group identified
    by \code{TRUE} values in \code{exp}. The length of \code{exp} and
    \code{obs} must be identical.}

  \item{\dots}{Additional arguments, available to methods.}
}

\value{

  A \code{data.frame} with columns

  \item{TruePositive}{Cummulative probability of correct assignment.}

  \item{FalsePositive}{Cummulative probability of incorrect assignment.}

}

\author{Martin Morgan \url{mailto:mtmorgan@fhcrc.org}}

\examples{
library(lattice)

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

## count data used for cross-validation, and cross-validation
count <- csubset(c("Lean", "Obese"), count, pheno)
data(bestgrp)

## true, false positives from single-group classifier
bst <- roc(pheno[rownames(count)] == "Obese",
           predict(bestgrp, count)[,"Obese"])
head(bst)

## lattice plot
xyplot(TruePostive ~ FalsePositive, bst, type="l",
       xlab="False Positive", ylab="True Positive")
}

\keyword{stats}
