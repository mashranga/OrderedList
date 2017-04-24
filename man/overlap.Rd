\name{overlap}
\alias{overlap}
\title{Count Elements in Overlap between two Lists}
\description{
  For each rank up to a given limit, count the number of elements in
  the overlap between two lists.
}
\usage{
overlap(x1, x2, n)
}
\arguments{
  \item{x1,x2}{ordered lists}
  \item{n}{the largest rank to be considered}
}
\value{
  Returns a vector of integers. The i-th element gives 
  the number of common elements in the first i positions 
  of both lists.
}
\author{Claudio Lottaz}
\keyword{internal}

