\name{sage.library}
\alias{read.sage.library}
\alias{write.sage.library}
\alias{plot.sage.library}
\alias{summary.sage.library}
\alias{print.sage.library}
\title{Class sage.library}
\description{
  The SAGE library class contains all the data and annotation for a SAGE
  library. It can contain two data.frames. 
}
\usage{
read.sage.library(file)
write.sage.library(x, file=paste(x\$libname, "sage", sep="."),
                   what="complete") 
}
\arguments{
  \item{x}{A sage library object}
  \item{file}{File name to read or write to}
  \item{what}{"complete", read complete librarary tags and sequences;
    "tags", read only tags and counts}
}
\details{
  SAGE library objects consists of one or two data.frames. The
  data.frame "tags" contains all the unique tags in the library and its
  counts. The data.frame "seqs" contains all the individual tag
  sequences and associated quality values. \code{read.sage.library} and
  \code{write.sage.library} are utility functions to read and write SAGE
  libraries.
  
}
\references{\url{http://tagcalling.mbgproject.org}}
\author{Tim Beissbarth}
\seealso{\code{\link[sagenhaft:extract.lib]{extract.lib}}}
\examples{
library(sagenhaft)
E15postHFI <- read.sage.library(system.file("data/E15postHFI.sage",
                                package="sagenhaft")) 
E15postHFI
}
\keyword{IO}
\concept{SAGE}

