\name{kendall.u}
\alias{kendall.u}
\alias{print.kendall.u}
\title{Kendall's Coefficient of Agreement}
\description{
  Kendall's u coefficient of agreement between judges.
}
\usage{kendall.u(M, correct = TRUE)}
\arguments{
  \item{M}{a square matrix or a data frame consisting of absolute choice
    frequencies; row stimuli are chosen over column stimuli}
  \item{correct}{logical, if \code{TRUE} (default) a continuity correction is
    applied when computing the test statistic (by subtracting one from the sum
    of agreeing pairs)}
}
\details{
  Kendall's u (Kendall and Babington Smith, 1940) takes on values between
  \code{min.u} (minimum agreement) and 1 (maximum agreement).
  The minimum \code{min.u} equals \eqn{-1/(m - 1)}, if \eqn{m} is even,
  and \eqn{-1/m}, if \eqn{m} is odd, where \eqn{m} is the number of subjects
  (judges).

  The null hypothesis in the chi-square test is that the agreement between
  judges is by chance.

  It is assumed that there is an equal number of observations per pair
  and that each subject judges each pair only once.
}
\value{
  \item{u}{Kendall's u coefficient of agreement}
  \item{min.u}{the minimum value for u}
  \item{chi2}{the chi-square statistic for a test that the agreement is
    by chance}
  \item{df}{the degrees of freedom}
  \item{p.value}{the p-value of the test}
}
\references{
  Kendall, M.G., & Babington Smith, B. (1940).
  On the method of paired comparisons.
  \emph{Biometrika}, \bold{31}, 324--345.
  \doi{10.1093/biomet/31.3-4.324}
}
\seealso{
  \code{\link{schoolsubjects}}, \code{\link{eba}}, \code{\link{strans}},
  \code{\link{circular}}.
}
\examples{
data(schoolsubjects)
lapply(schoolsubjects, kendall.u)  # better-than-chance agreement
}
\keyword{models}
