\name{simulate.eba}
\alias{simulate.eba}
\title{Simulate Responses from EBA Models}
\description{
  Simulates responses from the distribution corresponding to a fitted
  \code{eba} model object.
}
\usage{
\method{simulate}{eba}(object, nsim, seed, pool = TRUE, \dots)
}
\arguments{
  \item{object}{an object of class \code{eba}, typically the result of a
    call to \code{\link{eba}}}
  \item{nsim}{currently not used}
  \item{seed}{currently not used}
  \item{pool}{logical, if TRUE (default), pooled responses (summed across
    respondents) are returned}
  \item{\dots}{further arguments passed to or from other methods.
    None are used in this method.}
}
\details{
  Responses are simulated by \code{\link{rbinom}} with sizes taken from the
  original sample and probabilities computed from the model object.
}
\value{
  A paired-comparison matrix of (pooled) responses.
}
\seealso{
  \code{\link{eba}}, \code{\link{boot}}.
}
\examples{
data(celebrities)                     # absolute choice frequencies
A <- list(c(1,10), c(2,10), c(3,10),
          c(4,11), c(5,11), c(6,11),
          c(7,12), c(8,12), c(9,12))  # the structure of aspects
eba1 <- eba(celebrities, A)           # fit elimination-by-aspects model

## Perform parametric bootstrap
LR.stat <- replicate(200, eba(simulate(eba1), A)$good[1])

hist(LR.stat, col="lightgrey", border="white", freq=FALSE, breaks=20,
  xlim=c(0, 60), main="Parametric bootstrap")
curve(dchisq(x, df=eba1$good[2]), add=TRUE)
abline(v=eba1$good[1], lty=2)
}
\keyword{models}
