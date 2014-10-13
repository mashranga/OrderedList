\name{preparePermutations}
\alias{preparePermutations}
\title{Prepare Permutation and Subsetting Matrices}
\description{
  For a dataset specified with matrix and class labels, draws permutations
  and subsets.
}
\usage{
preparePermutations(ids, paired, B, sample.ratio = 0.8)
}
\arguments{
  \item{ids}{class labels}
  \item{paired}{logical, whether samples in classes are paired.}
  \item{B}{number of permutations and subsets to be drawn.}
  \item{sample.ratio}{how many of the samples in a class are to be subsampled.}
}
\value{
  Returns a list with the following items:
  \item{yperm}{the matrix of permutations.}
  \item{ysubs}{the matrix of subsamplings.}
}
\author{Claudio Lottaz}
\keyword{internal}
