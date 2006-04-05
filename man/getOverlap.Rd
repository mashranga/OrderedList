\name{getOverlap}
\alias{getOverlap}
\alias{print.listComparisonOverlap}
\alias{plot.listComparisonOverlap}
\title{ Extracting the Intersecting IDs From a listComparison Object }
\description{
This function extracts the intersecting set of list identifiers from an object of class \code{listComparison} as output of function \code{compareLists}. The user has to specify the maximum rank to be considered to receive the intersecting set up to this rank.
}
\usage{
getOverlap(x, max.rank = 100, percent = 0.95)
}
\arguments{
  \item{x}{ An object of class \code{listComparison}. }
  \item{max.rank}{ The maximum rank to be considered. }
  \item{percent}{ The final list of overlapping genes consists of those probes that contribute a certain percentage to the overall similarity score. Default is \code{percent=0.95}. To get the full list of genes, set \code{percent=1}. }
}
\details{
Function \code{compareLists} returns a list comparison for several choices of alpha. The number of genes to be taken into account differs dependent on alpha. One might now want to fix the number of genes and receive the resulting set of intersecting list identifiers. To this end, the user chooses a maximum rank to be considered from the values in column 'Genes' of the \code{listComparison} object. The direction (original or reversed) will internally be set to the direction yielding the higher similarity score.

If \code{two.sided} was \code{TRUE}, the first \code{max.rank} IDs on top of the lists and the first \code{max.rank} identifiers at the end of the lists are considered. If \code{two.sided} was \code{FALSE}, only the \code{max.rank} top identifiers are evaluated for overlap.
}
\value{
  An object of class \code{listComparisonOverlap} is returned. It contains the following
  list elements:
  \item{n}{the length of the lists.}
  \item{call}{the parameters of the input object.}
  \item{nn}{the input max.rank.}
  \item{score}{the observed similarity score.}
  \item{pvalue}{p-values for the observed score.}
  \item{overlaps}{number of overlapping identifiers per rank.}
  \item{randomScores}{random scores for given parameters.}
  \item{direction }{numerical value. Returns '1' if the similarity score is higher for the originally ordered lists and '-1' if the score is higher for the comparison of one original to one reversed list.}
  \item{intersect }{Vector with the sorted overlapping list identifiers, which contribute \code{percent} to the overall similarity score.}

  There are print and plot methods for \code{listComparisonOverlap} objects. The plot
  method takes a parameter \code{which} to specify whether "overlap" or
  "density" is to be drawn.
}
\references{
Yang X, Bentink S, Scheid S, and Spang R (2006): Similarities of ordered gene lists, to appear in \emph{Journal of Bioinformatics and Computational Biology}.
}
\author{ Claudio Lottaz, Stefanie Scheid }
\seealso{ \code{\link{OrderedList}}, \code{\link{compareLists}} }
\examples{
### Compare two artificial lists with some overlap
data(OL.data)
list1 <- as.character(OL.data$map$prostate)
list2 <- c(sample(list1[1:500]),sample(list1[501:1000]))
x <- compareLists(list1,list2)
x
getOverlap(x)
}
\keyword{datagen}
