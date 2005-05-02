\name{pcX}
\alias{pcX}
\title{Paired-Comparison Design Matrix}
\description{
  Computes a paired-comparison design matrix.
}
\usage{pcX(nstimuli)}
\arguments{
  \item{nstimuli}{number of stimuli in the paired-comparison design}
}
\details{
  The design matrix can be used when fitting a Bradley-Terry-Luce (BTL)
  model by means of \code{\link{glm}}.
}
\value{
  A matrix having \code{(nstimuli-1)*nstimuli/2} rows and
  \code{nstimuli} columns.
}
\seealso{\code{\link{OptiPt}}, \code{\link{glm}}.}
\examples{
data(celebrities)  # absolute choice frequencies
btl <- OptiPt(celebrities)
summary(btl)

btl.glm <- glm(cbind(btl$y1, btl$y0) ~ pcX(9)[,-1] - 1, binomial)
summary(btl.glm)  # BTL model does not fit these data!
}
\keyword{models}
