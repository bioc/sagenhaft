\name{error.correction}
\alias{estimate.errors.mean}
\alias{compute.sequence.neighbors}
\alias{em.estimate.error.given}

\title{Estimate sequencing errors and compute corrected counts}
\usage{
estimate.errors.mean(lib)
compute.sequence.neighbors(tags, taglength=10, quality.scores=NULL,
                           output="character") 
em.estimate.error.given(lib, maxstep=50, ...)
}
\description{
  These functions are used to compute sequencing error correction in a
  library. They are automatically called when extracting tags from
  sequences and therefore usually do not have to be called directly.
  
}
\arguments{
  \item{lib}{A sage library object}
  \item{tags}{A character vector or numeric vector containing tags}
  \item{taglength}{length of tag}
  \item{quality.scores}{A matrix containing base quality scores as -10
    log10 Pe}
  \item{maxstep}{iterations of EM algorithm}
  \item{output}{Output type for \code{compute.sequence.neighbors},
    either character or numeric}
  \item{\dots}{Other arguments ignored.}
}
\references{\url{http://tagcalling.mbgproject.org}}
\author{Tim Beissbarth}

\seealso{\code{\link[sagenhaft:extract.lib]{extract.lib}},
  \code{\link[sagenhaft:sage.library]{sage.library}}} 

\examples{
library(sagenhaft)
B6Hypo <-read.sage.library(system.file("data/B6HypothalHFI.sage",
                           package="sagenhaft")) 
E15post <- read.sage.library(system.file("data/E15postHFI.sage",
                             package="sagenhaft")) 
testlib <- combine.libs(B6Hypo, E15post)
testlib <- estimate.errors.mean(testlib)
testlib <- em.estimate.error.given(testlib)
tagneighbors <- compute.sequence.neighbors(testlib$seqs[,"seq"], 10,
                          testlib$seqs[,paste("q", 1:10, sep="")]) 
}
\keyword{error}
