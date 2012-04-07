\name{DMNGroup-class}
\Rdversion{1.1}
\docType{class}
\alias{DMNGroup-class}

\title{Class \code{"DMNGroup"}}

\description{Result from fitting a Dirichlet-Multinomial generative classifier.}

\section{Objects from the Class}{
  Objects can be created by calls to \code{\link{dmngroup}}.
}

\section{Slots}{

  All slots in this class are inheritted from \code{SimpleList}; see
  \sQuote{Methods}, below, for information on how to manipulate this
  object.

}

\section{Extends}{
  Class \code{"\linkS4class{SimpleList}"}, directly.
  Class \code{"\linkS4class{List}"}, by class "SimpleList", distance 2.
  Class \code{"\linkS4class{Vector}"}, by class "SimpleList", distance 3.
  Class \code{"\linkS4class{Annotated}"}, by class "SimpleList", distance 4.
}

\section{Methods}{

  See the \code{\link{mixture}} help page for functions that operate on
  \code{DMNGroup} and \code{DMN}.

  \code{DMNGroup} can be manipulated as a list; see
  \code{\linkS4class{SimpleList} } for a description of typical
  list-like functions.

}
\author{
  Martin Morgan \url{mailto:mtmorgan@fhcrc.org}
}

\seealso{
  \code{\link{mixture}}, \code{\linkS4class{DMN}},
  \code{\linkS4class{SimpleList}}. 
}
\examples{
data(bestgrp)
bestgrp
bestgrp[[1]]
}
\keyword{classes}
