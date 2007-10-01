\name{prepareData}
\alias{prepareData}
\title{ Combining Two Studies into an Expression Set }
\description{
The function prepares a collection of two expression sets (\code{ExpressionSet}) and/or Affy batches (\code{AffyBatch}) to be passed on to the main function \code{\link{OrderedList}}. For each data set, one has to specify the variable in the corresponding phenodata from which the grouping into two distinct classes is done. The data sets are then merged into one \code{ExpressionSet} together with the rearranged phenodata. If the studies were done on different platforms but a subset of genes can be mapped from one chip to the other, this information can be provided via the \code{mapping} argument.

Please note that both data sets have to be \emph{pre-processed} beforehand, either together or independent of each other. In addition, the gene expression values have to be on an \emph{additive scale}, that is logarithmic or log-like scale.
}
\usage{
prepareData(eset1, eset2, mapping = NULL)
}
\arguments{
  \item{eset1}{ The main inputs are the distinct studies. Each study is stored in a named list, which has five elements: \code{data}, \code{name}, \code{var}, \code{out} and \code{paired}, see details below. }
  \item{eset2}{ Same as \code{eset2} for the second data set. }
  \item{mapping}{ Data frame containing one named vector for each study. The vectors are comprised of probe IDs that fit to the rownames of the corresponding expression set. For each study, the IDs are ordered identically. For example, the \eqn{k}th row of \code{mapping} provides the label of the \eqn{k}th gene in each single study. If all studies were done on the same chip, no mapping is needed (default). }
}
\details{
Each study has to be stored in a list with five elements:
\tabular{ll}{
\code{data} \tab Object of class \code{ExpressionSet} or \code{AffyBatch}. \cr
\code{name} \tab Character string with comparison label. \cr
\code{var} \tab Character string with phenodata variable. Based on this variable, the samples for the two-sample testing will be extracted. \cr
\code{out} \tab Vector of two character strings with the levels of \code{var} that define the two clinical classes. The order of the two levels must be identical for all studies. Ideally, the first entry corresponds to the \emph{bad} and the second one to the \emph{good} outcome level. \cr
\code{paired} \tab Logical - \code{TRUE} if samples are paired (e.g. two measurements per patients) or \code{FALSE} if all samples are independent of each other. If data are paired, the paired samples need to be in (whatever) successive order. Thus, the first sample of one condition must match to the first sample of the second condition and so on.
}
}
\value{
An object of class \code{ExpressionSet} containing the joint data sets with appropriate phenodata.
}
\author{ Stefanie Scheid }
\references{
Yang X, Bentink S, Scheid S, and Spang R (2006): Similarities of ordered gene lists, to appear in \emph{Journal of Bioinformatics and Computational Biology}.
}
\seealso{ \code{\link{OL.data}}, \code{\link{OrderedList}} }
\examples{
data(OL.data)

### 'map' contains the appropriate mapping between 'breast' and 'prostate' IDs.
### Let's first concatenate two studies.
A <- prepareData(
                 list(data=OL.data$prostate,name="prostate",var="outcome",out=c("Rec","NRec"),paired=FALSE),
                 list(data=OL.data$breast,name="breast",var="Risk",out=c("high","low"),paired=FALSE),
                 mapping=OL.data$map
                 )

### We might want to examine the first 100 probes only.
B <- prepareData(
                 list(data=OL.data$prostate,name="prostate",var="outcome",out=c("Rec","NRec"),paired=FALSE),
                 list(data=OL.data$breast,name="breast",var="Risk",out=c("high","low"),paired=FALSE),
                 mapping=OL.data$map[1:100,]
                 )
}
\keyword{ datagen }
