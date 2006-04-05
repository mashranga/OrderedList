\name{OL.data}
\alias{OL.data}
\docType{data}
\title{ Gene Expression and Clinical Information of Two Cancer Studies }
\description{
The data contains a list with three elements: \code{breast}, \code{prostate} and \code{map}. The first two are expression sets of class \code{exprSet} taken from the breast cancer study of Huang et al. (2003) and the prostate cancer study of Singh et al. (2002). Both data sets were preprocessed as described in Yang et al. (2006). The data sets serve as illustration for function \code{\link{prepareData}}. Hence the sets contain only a random subsample of the original probes. We further removed unneeded samples from both studies.

The labels of the \code{breast} expression set were extended with 'B' to create two data sets where the probe IDs differ but can be mapped onto each other. The mapping is stored in the data frame \code{map}, which consists of the two probe ID vectors.
}
\usage{data(OL.data)}
\references{
Huang E, Cheng S, Dressman H, Pittman J, Tsou M, Horng C, Bild A, Iversen E, Liao M, Chen C, West M, Nevins J, and Huang A (2003): Gene expression predictors of breast cancer outcomes, \emph{Lancet} \bold{361}, 1590--1596.

Singh D, Febbo PG, Ross K, Jackson DG, Manola J, Ladd C, Tamayo P, Renshaw AA, D'Amico AV, Richie JP, Lander E, Loda M, Kantoff PW, Golub TR, and Sellers WR (2002): Gene expression correlates of clinical prostate cancer behavior, \emph{Cancer Cell} \bold{1}, 203--209.

Yang X, Bentink S, Scheid S, and Spang R (2006): Similarities of ordered gene lists, to appear in \emph{Journal of Bioinformatics and Computational Biology}.
}
\seealso{ \code{\link{OL.result}} }
\keyword{datasets}



