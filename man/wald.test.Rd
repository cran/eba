\name{wald.test}
\alias{wald.test}
\alias{print.wald.test}
\title{Testing Linear Hypotheses}
\description{
  Tests linear hypotheses of the form \eqn{Cp = 0} in EBA models using
  the Wald test.
}
\usage{wald.test(object, C, u.scale = TRUE)}
\arguments{
  \item{object}{an object of class \code{eba}, typically the result of a
    call to \code{\link{eba}}}
  \item{C}{a matrix of contrasts, specifying the linear hypotheses}
  \item{u.scale}{logical, if TRUE the test is performed on the u-scale, if
                 FALSE the test is performed on the EBA parameters directly}
}
\details{
  The Wald test statistic,
  \deqn{W = (Cp)' [C cov(p) C']^{-1} (Cp),}
  is approximately chi-square distributed with \eqn{rk(C)} degrees of
  freedom.

  \code{C} is usually of full rank and must have as many columns as there
  are parameters in \code{p}.
}
\value{
  \item{C}{the matrix of contrasts, specifying the linear hypotheses}
  \item{W}{the Wald test statistic}
  \item{df}{the degrees of freedom (\eqn{rk(C)})}
  \item{pval}{the p-value of the test}
}
\seealso{
  \code{\link{eba}}, \code{\link{group.test}}, \code{\link{cov.u}}.
}
\examples{
data(celebrities)  # absolute choice frequencies
A <- list(c(1,10), c(2,10), c(3,10), c(4,11), c(5,11), c(6,11),
          c(7,12), c(8,12), c(9,12))  # the structure of aspects
eba1 <- eba(celebrities, A)  # fit a preference tree

## Test whether JU, CY, and AJF have equal preference scale values
C1 <- matrix(c(0,0,0,1,-1,0,0,0,0,
               0,0,0,1,0,-1,0,0,0), 2, 9, TRUE)
wald.test(eba1, C1)

## Test whether the three branch parameters are different
C2 <- matrix(c(0,0,0,0,0,0,0,0,0,1,-1,0,
               0,0,0,0,0,0,0,0,0,1,0,-1), 2, 12, TRUE)
wald.test(eba1, C2, u.scale = FALSE)
}
\keyword{models}
