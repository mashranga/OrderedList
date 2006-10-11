\name{scoreOrderComparison}
\alias{scoreOrderComparison}
\alias{scoreOrderComparisonBoth}
\title{Score the Comparison of two Gene Rankings}
\description{
  Compute weighted similarity score for gene rankings determined via
  the chosen test statistics.
}
\usage{
scoreOrderComparison(exprs1, labels1, paired1, 
                     exprs2, labels2, paired2, 
                     test.method = test.z, nn, bases, two.sided, empirical)
scoreOrderComparisonBoth(exprs1, labels1, paired1, 
                         exprs2, labels2, paired2, 
                         test.method = test.z, nn, bases, two.sided, empirical)
}
\arguments{
  \item{exprs1,exprs2}{gene expression matrices.}
  \item{labels1,labels2}{class labels, one label per column in matrices.}
  \item{paired1,paired2}{logical, whether samples are paired in classes.}
  \item{test.method}{a function computing one test statistics per row and taking
    a matrix, a label vector and a logical for pairing as parameters. 
    Valid examples are \code{test.fc}, \code{test.t} and \code{test.z}.}
  \item{nn}{a vector of rank limits. The score is computed taking into
    account ranks up to these limits only. One limit per entry in \code{bases}.}
  \item{bases}{a vector of bases used in weighted scores, is equal to exp(-alpha).
    The function can compute scores for several regularization parameters in one go.}
  \item{two.sided}{if \code{TRUE} both ends of the lists are taken into account, only
    top ranks are considered otherwise.}
}
\details{
  \code{scoreOrderComparison} computes scores only for the direct
  comparison. \code{scoreOrderComparisonBoth} in addition computes
  scores for reversed orderes, i.e., one of the rankings is reversed.
}
\value{
  For each entry in \code{bases}, thus for each regularization parameter alpha,
  one score is returned in an array.
}
\author{Claudio Lottaz}
\keyword{internal}
