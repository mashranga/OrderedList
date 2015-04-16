\name{shuffledRandomScores}
\alias{shuffledRandomScores}
\alias{plot.shuffledRandomScores}
\alias{print.shuffledRandomScores}
\title{Generates Null-Distribution for List-Overlap-Scores}
\description{
  A null-distribution for list-overlap scores is generated via 
  simulation. Scores are computed for random permutations. 
}
\usage{
shuffledRandomScores(n, nn, bases, B = 1000, two.sided=TRUE)
}
\arguments{
  \item{n}{the length of the lists.}
  \item{nn}{the maximal ranks to be considered, one entry per weighting parameter alpha.}
  \item{bases}{\code{exp(-alpha)} for each weighting parameter alpha.}
  \item{B}{the number of permutations to be drawn.}
  \item{two.sided}{if \code{TRUE} both ends of the lists are taken into account, only
    top ranks are considered otherwise.}
}
\value{
  Returns an object of type "shuffledRandomScores. Its only data is 
  a matrix of random scores. One column per alpha and one
  row per permutation is generated.

  There are print and plot methods for "shuffledRandomScores objects.
}
\author{Claudio Lottaz}
\seealso{\code{\link{compareLists}}}
\keyword{internal}

