\name{sage.library.comparison}
\alias{read.sage.library.comparison}
\alias{write.sage.library.comparison}
\alias{plot.sage.library.comparison}
\alias{summary.sage.library.comparison}
\alias{print.sage.library.comparison}
\alias{compare.lib.pair}
\title{Class sage.library.comparison}
\description{
  Class for storing the data of a pairwise comparison between two SAGE
  libraries.

}
\usage{
read.sage.library.comparison(file)
write.sage.library.comparison(x, file=paste(lib\$libname, "sage", sep="."))
compare.lib.pair(lib1, lib2) 
}
\arguments{
  \item{x,lib1,lib2}{A sage library object}
  \item{file}{File name to read or write to}
}
\details{
  
  SAGE library comparison objects consists of one data.frames. It stores
  a A and an M value which are the log2 average expression and log2
  ratio, respectively. It also has a column for the resulting p.values
  from \code{\link[sagenhaft:sage.test]{sage.test}}.
  \code{read.sage.library.comparison} and
  \code{write.sage.library.comparison} are utility functions to read and
  write SAGE library comparisons. \code{compare.lib.pair} can be used to
  generate SAGE library comparisons.
  
}
\references{\url{http://tagcalling.mbgproject.org}}
\author{Tim Beissbarth}
\seealso{
\code{\link[sagenhaft:sage.test]{sage.test}}
}
\examples{
library(sagenhaft)
B6Hypo <- read.sage.library(system.file("data/B6HypothalHFI.sage",
                            package="sagenhaft")) 
E15post <- read.sage.library(system.file("data/E15postHFI.sage",
                             package="sagenhaft")) 
libcomp <- compare.lib.pair(B6Hypo, E15post)
plot(libcomp)
libcomp
}
\keyword{IO}
\concept{SAGE}

