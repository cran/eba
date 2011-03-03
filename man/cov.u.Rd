\name{cov.u}
\alias{cov.u}
\title{Covariance Matrix of the EBA Utility Scale}
\description{
  Computes the (normalized) covariance matrix of the utility scale from the
  covariance matrix of elimination-by-aspects (EBA) model parameters.
}
\usage{cov.u(object, norm = "sum")}
\arguments{
  \item{object}{an object of class \code{eba}, typically the result of a
    call to \code{\link{eba}}}
  \item{norm}{either \code{sum} (default), a number from 1 to number of
    stimuli, or \code{NULL}; see \code{\link{uscale}} for details}
}
\details{
  The additivity rule for covariances
  \eqn{cov(x + y, z) = cov(x, z) + cov(y, z)} is used for the computations.

  If \code{norm} is not \code{NULL}, the unnormalized covariance matrix is
  transformed using \eqn{a^2 cov(u)}, where the constant \eqn{a} results from
  the type of normalization applied.
}
\value{
  The (normalized) covariance matrix of the utility scale.
}
\seealso{
  \code{\link{uscale}}, \code{\link{eba}}, \code{\link{wald.test}}.
}
\keyword{models}
