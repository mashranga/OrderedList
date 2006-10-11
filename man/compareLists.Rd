\name{compareLists}
\alias{compareLists}
\alias{print.listComparison}
\alias{plot.listComparison}
\title{Compare Ordered Lists with Weighted Overlap Score}
\description{
  The two orderings recieved as parameters are compared using the 
  weighted overlap score and compared with a random distribution 
  of that score (yielding an empirical p-value). 
}
\usage{
compareLists(ID.List1, ID.List2, mapping = NULL, 
             two.sided=TRUE, B = 1000, alphas = NULL, 
             invar.q = 0.5, min.weight = 1e-5)
}
\arguments{
  \item{ID.List1}{first ordered list of identifiers to be compared.}
  \item{ID.List2}{second ordered list to be compared, must have the same 
    length as \code{ID.List1}.}
  \item{mapping}{maps identifiers between the two lists. This is a matrix with 
    two columns. All items in \code{ID.List1} must match to exactly one entry of 
    column 1 of the mapping, each element in \code{ID.List2} must match exactly one 
    element in column 2 of the mapping. If mapping is \code{NULL}, the two lists 
    are expected to contain the same identifiers and there must be a one-to-one
    relationship between the two.}
  \item{two.sided}{whether the score is to be computed considering both ends of the list,
    or just the top members.}
  \item{B}{the number of permutations used to estimate empirical p-values.}
  \item{alphas}{a set of alpha candidates to be evaluated. If set to
    \code{NULL}, alphas are
    determined such that reasonable maximal ranks to be considered result.}
  \item{invar.q}{quantile of genes expected to be invariant. These are
    not used
    during shuffling, since they are expected to stay away from the ends of the 
    lists, even when the data is perturbed to generate the NULL
    distribution. The default of 0.5 is reasonnable for whole-genome gene
    expression analysis, but must be reconsidered when the compared lists
    are deduced from other sources.}
  \item{min.weight}{the minimal weight to be considered. }
}
\details{
  The two lists received as arguments are matched against each other
  according to the given mapping. The comparison is performed from both
  ends by default. Permutations of lists are used to generate random scores
  and compute empirical p-values. The evaluation is also performed for the
  case the lists should be reversed. From the resulting output, the set of
  overlapping list identifiers can be extracted using function \code{getOverlap}.
}
\value{
  An object of class \code{listComparison} is returned. It contains the following
  list elements:
  \item{n}{the length of the lists}
  \item{call}{the input parameters}
  \item{nn}{the maximal number of genes corresponding to the alphas and the minimal weight}
  \item{scores}{scores for the straight list comparisons}
  \item{revScores}{scores for the reversed list comparison}
  \item{pvalues}{p-values for the straight list comparison}
  \item{revPvalues}{p-values for the reversed list comparison}
  \item{overlap}{number of overlapping identifiers per rank in straight comparison}
  \item{revOverlap}{number of overlapping identifiers per rank in reversed comparison}
  \item{randomScores}{random scores per weighting parameter}
  \item{ID.List1}{same as input ID.List1}
  \item{ID.List2}{same as input ID.List2}

  There are print and plot methods for \code{listComparison} objects. The plot
  method takes a parameter \code{which} to specify whether "overlap" or
  "density" is to be drawn.
}
\references{
Yang X, Bentink S, Scheid S, and Spang R (2006): Similarities of ordered gene lists, to appear in \emph{Journal of Bioinformatics and Computational Biology}.
}
\author{Claudio Lottaz, Stefanie Scheid}
\seealso{ \code{\link{OrderedList}}, \code{\link{getOverlap}} }
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

