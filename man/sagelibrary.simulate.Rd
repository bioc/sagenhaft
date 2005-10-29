\name{sagelibrary.simulate}
\alias{sagelibrary.simulate}
\title{ Simulate SAGE libraries }
\description{
  Function to simulate SAGE libraries with sequencing errors.
}
\usage{
sagelibrary.simulate(taglength = 4, lambda = 1000, mean.error = 0.01,
                 error.sd = 1, withintagerror.sd = 0.2,
                 ngenes = min(4^taglength, 1e+05), base.lib = NULL,
                 libseed = -1, ...)
}
\arguments{
  \item{taglength}{ Tag length for library. }
  \item{lambda}{ Aproximate size of library. }
  \item{mean.error}{ Mean amount of sequencing errors. }
  \item{error.sd}{ Standard deviation for sequencing errors. }
  \item{withintagerror.sd}{ Standard deviation for sequencing
    errors within tags. }
  \item{ngenes}{ Number of genes to generate tags from. }
  \item{base.lib}{ Simulate library based on tags in other lib and
    create variations. }
  \item{libseed}{ Seed for random number generator. }
  \item{\dots}{ Arguments passed to em.estimate. }
}
\details{
  We set the number of possible transcripts and assign a random SAGE tag
  to each of them out of all 4\^taglength possible SAGE tags. For each
  SAGE tag a random proportion p within the library is generated from a
  log-normal distribution, and the proportions are then adjusted to have
  a sum of 1. The true counts of a tag are simulated by sampling from
  Poisson distributions with parameters p lambda, where p is the
  proportion of the tag in the library and lambda is a parameter for
  setting the size of the library. The simulation of the sequencing
  errors is done on each individual occurrence of a tag sequence. For
  each tag sequence a mean sequencing quality value is generated from a
  log-normal distribution. The individual quality values for each base
  are then generated from log-normal distributions with means equal to
  the simulated sequencing quality values for the tag sequences. We have
  noticed that with experimentally generated data the within tag
  sequence variation of sequencing quality values is usually about 1/5
  of the between tag sequence variation. From each true tag sequence one
  observed tag sequence is generated using the simulated quality values
  of the true sequence as the multinomial probabilities, i.e. replacing
  each base with either one of the 3 other bases with the probability
  specified by the sequencing quality value of that base. The counts of
  these generated tags are then summed to represent the observed
  tags. When generating several simulated libraries for comparisons, we
  use the same proportions of the genes for all libraries, replacing up
  to 1/3 of the proportions by proportions with a known differential
  factor.
}
\references{\url{http://tagcalling.mbgproject.org}}
\author{Tim Beissbarth}

\seealso{\code{\link[sagenhaft:sage.library]{sage.library}},
  \code{\link[sagenhaft:error.correction]{error.correction}}}
\examples{
library(sagenhaft)
testlib1 <- sagelibrary.simulate(taglength=10, lambda=10000,
                             mean.error=0.01)
testlib2 <- sagelibrary.simulate(taglength=10, lambda=20000,
                             mean.error=0.02, base.lib=testlib1)
testlib3 <- sagelibrary.simulate(taglength=10, lambda=10000,
                             mean.error=0.01, libseed=testlib1$seed)
}
\keyword{ misc }

