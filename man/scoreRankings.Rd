\name{scoreRankings}
\alias{scoreRankings}
\title{Score the Comparison of two Rankings}
\description{
  Two rankings are accepted as input in the form of corresponding ranks in two lists. The weighted overlap
  score is computed efficiently without explicitly computing overlaps.
}
\usage{
scoreRankings(r1, r2, nn, bases, two.sided=TRUE)
}
\arguments{
  \item{r1}{integer, ranks in the first list.}
  \item{r2}{integer, ranks in the second list.
    \code{r1} and \code{r2} must have the same length.}
  \item{nn}{for each alpha to be used as weighting parameter,
    this array of integers contains the maximal rank for which
    overlaps are considered.}
  \item{bases}{for each alpha to be used as weighting parameter,
    this double array contains exp(-alpha).}
  \item{two.sided}{if \code{TRUE} both ends of the lists are taken into account, only
    top ranks are considered otherwise.}
}
\details{
  The score to be computed is defined as the sum over the first ranks
  in two lists. The summed up measure is the weighted overlap between the
  two lists:

  \deqn{score := sum_(R=1)^n exp(-alpha R) * overlap(L1[1:R],L2[1:R])}

  where \eqn{n} is the maximal rank to be considered and \eqn{L1/L2} denote the sorted
  lists to be compared. In this score, each gene contributes from the first
  rank where it is in the overlap up to \eqn{n}. For gene i:

  \deqn{score_i = sum_(R=max(r1[i],r2[i])^n exp(-alpha R))}

  where \eqn{r1/r2} are the ranks of genes in \eqn{L1/L2}. Since this is a finite geometric series, it can be used to speed up computation of our score:

  \deqn{score_i = (exp(-alpha min(r1[i],r2[i])) - exp(-alpha n))/(1 - exp(-alpha))}
  \deqn{score = sum_(i|r1[i]<n&r2[i]<n) score_i}

  Analogue computations are performed by \code{scoreRankings} for list 
  begins and list ends.
}
\value{
  An array of  doubles with one score per weighting parameter to 
  be considered is returned.
}
\author{Claudio Lottaz}
\seealso{\code{\link{shuffledRandomScores}},\code{\link{compareLists}}}
\keyword{internal}
