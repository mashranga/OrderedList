\name{check.test.args}
\alias{check.test.args}
\alias{test.fc}
\alias{test.t}
\alias{test.z}
\title{Helper Function to Compute Test Statistics}
\description{
  These function compute test statistics for all rows in 
  a matrix. 
}
\usage{
test.fc(m, cl, paired)
test.t(m, cl, paired)
test.z(m, cl, paired)
check.test.args(m, cl, paired)
}
\arguments{
  \item{m}{the matrix of numeric values. For each row of the matrix,
    one test statistic is computed.}
  \item{cl}{a vector of class labels. \code{cl} must have as many elements
    as there are columns in \code{m}.}
  \item{paired}{logical, is \code{TRUE}, if the data in the two classes are paired.}
}
\details{
  \code{check.test.arg} is used by the other methods to check their arguments.
  \code{test.t} and \code{test.z} interface with the C-code contained in the 
  twilight package in order to speed-up the computation.
}
\value{
  An array of the corresponding test statistic containing one element per 
  row of the input matrix.
}
\author{Claudio Lottaz}
\seealso{\code{\link[twilight]{twilight.teststat}}}
\keyword{internal}
