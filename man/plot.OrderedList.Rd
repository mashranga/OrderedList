\name{plot.OrderedList}
\alias{plot.OrderedList}
\title{ Plotting Function for OrderedList Objects }
\description{
The function generates three different plots, which can be selected via argument \code{which}. With default \code{which=NULL}, all three figures are plotted into one graphics device.
}
\usage{
\method{plot}{OrderedList}(x, which = NULL, no.title=FALSE, \dots)
}
\arguments{
  \item{x}{ Object of class \code{OrderedList}. }
  \item{which}{ Select one of the three figures described in the details section below. }
  \item{no.title}{logical, whether to skip plotting a title.}
  \item{\dots}{ Additional graphical arguments. }
}
\details{
\code{which} is one of 'pauc', 'scores' or 'overlap'. If \code{NULL}, all figures are produced in a row.

Option 'pauc' selects the plot of pAUC-scores, based on which the optimal \eqn{\alpha} is chosen. The pAUC-score measure the separability between the two distributions of observed and expected similarity scores. The similarity scores depend on \eqn{\alpha} and thus \eqn{\alpha} is chosen where the pAUC-scores are maximal. The optimal \eqn{\alpha} is marked by a vertical line.

Figure 'scores' shows kernel density estimates of the two score distributions underlying the pAUC-score for optimal \eqn{\alpha}. The red curve correspondence to simulated observed scores and the black curve to simulated expected scores. The vertical red line denotes the actually observed similarity score. The bottom rugs mark the simulated values. The two distributions got the highest pAUC-score of separability and thus provide the best signal-to-noise separation.

Finally, 'overlap' displays the numbers of overlapping genes in the two gene lists. The overlap size is drawn as a step function over the respective ranks. Top ranks correspond to up-regulated and bottom ranks to down-regulated genes. In addition, the expected overlap and 95\% confidence intervals derived from a hypergeometric distribution are plotted. If \code{empirical=TRUE} in \code{\link{OrderedList}} the confidence intervals were derived empirically from shuffling the data and computing the overlap under the null hypothesis.}
\value{
No value is returned.
}
\references{
Yang X, Bentink S, Scheid S, and Spang R (2006): Similarities of ordered gene lists, to appear in \emph{Journal of Bioinformatics and Computational Biology}.
}
\author{ Xinan Yang, Stefanie Scheid }
\seealso{ \code{\link{OrderedList}} }
\examples{
data(OL.result)
plot(OL.result)
}
\keyword{ hplot }
