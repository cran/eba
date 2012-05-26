\name{logLik.eba}
\alias{logLik.eba}
\alias{nobs.eba}
\title{Log-Likelihood of an eba Object}
\description{
  Returns the log-likelihood value of the (multi-attribute) probabilistic
  choice model represented by \code{object} evaluated at the
  estimated parameters.
}
\usage{
\method{logLik}{eba}(object, \dots)
}
\arguments{
 \item{object}{an object inheriting from class \code{eba}, representing
   a fitted elimination-by-aspects model.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\value{
  The log-likelihood of the model represented by
  \code{object} evaluated at the estimated parameters.
}
\seealso{
  \code{\link{eba}}, \code{\link{logLik.lm}}, \code{\link{AIC}}.
}
\examples{
data(heaviness)
btl1 <- eba(heaviness[,,order=1])
logLik(btl1)
AIC(btl1)
BIC(btl1)
}
\keyword{models}
