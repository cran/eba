\name{cov.u}
\alias{cov.u}
\title{Covariance of the u-scale}
\description{
  Computes the covariance matrix of the u-scale from the covariance
  matrix of the parameters.
}
\usage{cov.u(eba)}
\arguments{
  \item{eba}{an object of class \code{eba}}
}
\details{
  The additivity law of covariances \eqn{cov(x+y,z) = cov(x,z)+cov(y,z)}
  is used for the computations.
}
\value{
  the covariance matrix of the u-scale
}
\author{Florian Wickelmaier \email{wickelmaier@web.de}}
\seealso{
  \code{\link{OptiPt}}.
}
\keyword{models}
