\name{cov.u}
\alias{cov.u}
\title{Covariance of the u-scale}
\description{
  Computes the covariance matrix of the u-scale from the covariance
  matrix of the EBA parameters.
}
\usage{cov.u(eba)}
\arguments{
  \item{eba}{an object of class \code{eba}, typically the result of a
    call to \code{\link{OptiPt}}}
}
\details{
  The additivity law of covariances \eqn{cov(x+y,z) = cov(x,z)+cov(y,z)}
  is used for the computations.
}
\value{
  The covariance matrix of the u-scale.
}
\seealso{
  \code{\link{OptiPt}}, \code{\link{wald.test}}.
}
\keyword{models}
