\name{OrderedList}
\alias{OrderedList}
\title{ Detecting Similarities of Two Microarray Studies }
\description{
  Function \code{OrderedList} aims for the \emph{comparison of
  comparisons}: given two expression studies with one ranked (ordered)
  list of genes each, we might observe considerable overlap among the
  top-scoring genes. \code{OrderedList} quantifies this overlap by
  computing a weighted similarity score, where the top-ranking genes
  contribute more to the score than the genes further down the list. The
  final list of overlapping genes consists of those probes that
  contribute a certain percentage to the overall similarity score.
}
\usage{
OrderedList(eset, B = 1000, test = "z", beta = 1, percent = 0.95, verbose = TRUE, alpha=NULL, min.weight=1e-5)
}
\arguments{
  \item{eset}{ Expression set containing the two studies of interest. Use \code{\link{prepareData}} to generate \code{eset}. }
  \item{B}{ Number of internal sub-samples needed to optimize alpha. }
  \item{test}{ String, one of 'fc' (log ratio = log fold change), 't' (t-test with equal variances) or 'z' (t-test with regularized variances). The z-statistic is implemented as described in Efron et al. (2001). }
  \item{beta}{ Either 1 or 0.5. In a comparison where the class labels of the studies match, we set \code{beta=1}. For example, in each single study the first class relates to bad prognosis while the second class relates to good prognosis. If a matching is not possible, we set \code{beta=0.5}. For example, we compare a study with good/bad prognosis classes to a study, in which the classes are two types of cancer tissues. }
  \item{percent}{ The final list of overlapping genes consists of those probes that contribute a certain percentage to the overall similarity score. Default is \code{percent=0.95}. To get the full list of genes, set \code{percent=1}. }
  \item{verbose}{ Logical value for message printing. }
  \item{alpha}{A vector of weighting parameters. If set to NULL (the default),
    parameters are computed such that top 100 to the top 2500 ranks receive
    weights above \code{min.weight}.}
  \item{min.weight}{The minimal weight to be taken into account while computing scores.}
}
\details{
In short, the similarity measure is computed as follows: Based on two-sample test statistics like the t-test, genes within each study are ranked from most up-regulated down to most down-regulated. Thus we have one ordered list per study. Now for each rank going both from top (up-regulated end) and from bottom (down-regulated end) we count the number of overlapping genes. The total overlap \eqn{A_n} for rank \eqn{n} is defined as:
\deqn{A_n = O_n (G_1,G_2) + O_n(f(G_1),f(G_2))}
where \eqn{G_1} and \eqn{G_2} are the two ordered list, \eqn{f(G_1)} and \eqn{f(G_2)} are the two flipped lists with the down-regulated genes on top and \eqn{O_n} is the size of the overlap of its two arguments. A preliminary version of the weighted overlap over all ranks \eqn{n} is then given as:
\deqn{T_\alpha(G_1,G_2) = \sum_n \exp{-\alpha n} A_n.}
The final similarity score includes the case that we cannot match the classes in each study exactly and thus do not know whether up-regulation in one list corresponds to up- or down-regulation in the other list. Here parameter \eqn{\beta} comes into play:
\deqn{
S_\alpha(G_1,G_2) = \max{ \beta T_\alpha(G_1,G_2), (1-\beta) T_\alpha (G_1,f(G_2)) }.
}
Parameter \eqn{\beta} is set by the user but parameter \eqn{\alpha} has to be tuned in a simulation using sub-samples and permutations of the original class labels.
}
\value{
Returns an object of class \code{OrderedList}, which consists of a list with entries:
  \item{n}{Total number of genes.}
  \item{label }{The concatenated study labels as provided by \code{eset}.}
  \item{p }{The p-value specifying the significance of the similarity.}
  \item{intersect }{Vector with sorted probe IDs of the overlapping genes, which contribute \code{percent} to the overall similarity score.}
  \item{alpha }{The optimal regularization parameter alpha.}
  \item{direction }{Numerical value. Returns '1' if the similarity score is higher for the originally ordered lists and '-1' if the score is higher for the comparison of one original to one flipped list. Of special interest if \code{beta=0.5}.}
  \item{scores }{Matrix of observed test scores with genes in rows and studies in columns.}
  \item{sim.scores }{List with four elements with output of the resampling with optimal \code{alpha}. \code{SIM.observed}: The observed similarity sore. \code{SIM.alternative}: Vector of observed similarity scores simulated using sub-sampling within the distinct classes of each study. \code{SIM.random}: Vector of random similarity scores simulated by randomly permuting the class labels of each study. \code{subSample}: \code{TRUE} to indicate that sub-sampling was used.}
  \item{pauc }{Vector with pAUC-scores for each candidate of the regularization parameter \eqn{\alpha}. The maximal pAUC-score defines the optimal \eqn{\alpha}. See also \code{\link{plot.OrderedList}}.}
  \item{call }{List with some of the input parameters.}
}
\references{
Yang X, Bentink S, Scheid S, and Spang R (2006): Similarities of ordered gene lists, to appear in \emph{Journal of Bioinformatics and Computational Biology}.

Efron B, Tibshirani R, Storey JD, and Tusher V (2001): Empirical Bayes analysis of a microarray experiment, \emph{Journal of the American Statistical Society} \bold{96}, 1151--1160.
}
\author{ Xinan Yang, Claudio Lottaz, Stefanie Scheid }
\seealso{ \code{\link{prepareData}}, \code{\link{OL.data}}, \code{\link{OL.result}}, \code{\link{plot.OrderedList}}, \code{\link{print.OrderedList}}, \code{\link{compareLists}}}
\examples{
### Let's compare the two example studies.
### The first entries of 'out' both relate to bad prognosis.
### Hence the class labels match between the two studies
### and we can use 'OrderedList' with default 'beta=1'.
data(OL.data)
a <- prepareData(
                 list(data=OL.data$breast,name="breast",var="Risk",out=c("high","low"),paired=FALSE),
                 list(data=OL.data$prostate,name="prostate",var="outcome",out=c("Rec","NRec"),paired=FALSE),
		 mapping=OL.data$map
                 )
\dontrun{
OL.result <- OrderedList(a)
}

### The same comparison was done beforehand.
data(OL.result)
OL.result
plot(OL.result)
}
\keyword{ htest }
