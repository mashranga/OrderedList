\name{OL.result}
\alias{OL.result}
\docType{data}
\title{ Three Examples of Class 'OrderedList' }
\description{
The data set consists of an \code{OrderedList} object derived by applying function \code{\link{OrderedList}} on the expression sets in \code{\link{OL.data}}. The function calls are given in the example section below.
}
\usage{data(OL.result)}
\seealso{ \code{\link{OL.data}}, \code{\link{OrderedList}} }
\references{
Yang X, Bentink S, Scheid S, and Spang R (2006): Similarities of ordered gene lists, to appear in \emph{Journal of Bioinformatics and Computational Biology}.
}
\examples{
\dontrun{
a <- prepareData(
                 list(data=OL.data$breast,name="breast",var="Risk",out=c("high","low"),paired=FALSE),
                 list(data=OL.data$prostate,name="prostate",var="outcome",out=c("Rec","NRec"),paired=FALSE),
		 mapping=OL.data$map
                 )
OL.result <- OrderedList(a)
}
}
\keyword{datasets}



