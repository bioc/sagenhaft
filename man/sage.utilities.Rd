\name{sage.utilities}
\alias{tagnum2tagmatrix}
\alias{tagmatrix2tagnum}
\alias{tagnum2tagsequence}
\alias{tagsequence2tagnum}
\alias{revcomp}
\alias{create.matrix.csr}
\alias{table.sparse}
\alias{difference.scatter.plot}

\title{Utilities}
\description{
  Different utilities to use with SAGE data.
}

\usage{
tagnum2tagmatrix(tags, length)
tagmatrix2tagnum(tags, length=ncol(tags))
tagnum2tagsequence(tags, length)
tagsequence2tagnum(tags, length)
revcomp(seq)
}

\arguments{
  \item{tags}{integer or character vector giving SAGE tags.}
  \item{length}{Length of SAGE tags.}
  \item{seq}{Character vector or list of sequences.}
  \item{...}{SAGE library objects.}
}

\details{
  These functions are utility functions used in SAGE tag extraction,
  e.g. to convert SAGE tag sequences to numeric values, i.e. base 4 for
  efficient storage and handling, and to reverse complement sequences.
}

\author{Tim Beissbarth}

\examples{
library(sagenhaft)
tags <- c("aaa", "ttt", "ccc")
tagsnumeric <- tagsequence2tagnum(tags, 3)
tagsmatrix <- tagnum2tagmatrix(tagsnumeric, 3)
tags <- tagnum2tagsequence(tagmatrix2tagnum(tagsmatrix, 3), 3)
revcomp(tags)
}

\keyword{utilities}
