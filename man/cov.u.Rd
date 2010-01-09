\name{cov.u}
\alias{cov.u}
\title{Covariance Matrix of the Utility Scale}
\description{
  Computes the covariance matrix of the utility scale from the covariance
  matrix of EBA parameters.
}
\usage{cov.u(object)}
\arguments{
  \item{object}{an object of class \code{eba}, typically the result of a
    call to \code{\link{eba}}}
}
\details{
  The additivity rule for covariances \eqn{cov(x+y,z) = cov(x,z) + cov(y,z)}
  is used for the computations.
}
\value{
  The covariance matrix of the utility scale.
}
\seealso{
  \code{\link{eba}}, \code{\link{wald.test}}.
}
\keyword{models}
