\name{kendall.u}
\alias{kendall.u}
\alias{print.kendall.u}
\title{Kendall's Coefficient of Agreement}
\description{
  Kendall's u coefficient of agreement between judges.
}
\usage{kendall.u(M, cont.correct = FALSE)}
\arguments{
  \item{M}{a square matrix or a data frame consisting of absolute choice
    frequencies; row stimuli are chosen over column stimuli}
  \item{cont.correct}{logical, if TRUE a correction for continuity is
    applied (by deducting 1 from \code{chi2}); default is FALSE}
}
\details{
  Kendall's u takes values between \code{min.u} (when agreement is minimum)
  and 1 (when agreement is maximum).
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
  \item{pval}{the p-value of the test}
}
\references{
  Kendall, M.G., & Smith, B.B. (1940). On the method of paired comparisons. 
  \emph{Biometrika}, \bold{31}, 324--345.
}
\seealso{
  \code{\link{eba}}, \code{\link{strans}}.
}
\examples{
data(celebrities)
kendall.u(celebrities)  # moderate agreement
}
\keyword{models}
