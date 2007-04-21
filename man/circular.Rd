\name{circular}
\alias{circular}
\alias{print.circular}
\title{Circular triads (intransitive cycles)}
\description{
  Number of circular triads and indices of inconsistency.
}
\usage{circular(mat)}
\arguments{
  \item{mat}{a square matrix or a data frame consisting of (individual)
    binary choice data; row stimuli are chosen over column stimuli}
}
\details{
  Kendall's coefficient of consistency,
  \deqn{zeta = 1 - T/T_{max},}
  lies between one (perfect consistency) and zero,
  where \code{T} is the observed number of cicular triads,
  and the maximum possible number of circular triads is
  \eqn{T_{max} = n*(n^2 - 4)/24}, if \eqn{n} is even, and
  \eqn{T_{max} = n*(n^2 - 1)/24} else, and \eqn{n} is the
  number of stimuli or objects judged.

  The null hypothesis in the chi-square test is that the inconsistencies
  in the paired-comparison judgments are by chance. The chi-square
  approximation might be incorrect if \eqn{n < 8}.
}
\value{
  \item{T}{number of cicular triads}
  \item{T.max}{maximum possible number of circular triads}
  \item{zeta}{Kendall's coefficient of consistency}
  \item{chi2}{the chi-square statistic for a test that the consistency is
    by chance}
  \item{df}{the degrees of freedom}
  \item{pval}{the p-value of the test}
}
\references{
  David, H. (1988). \emph{The method of paired comparisons}. 
  London: Griffin.
}
\seealso{
  \code{\link{eba}}, \code{\link{strans}}, \code{\link{kendall.u}}.
}
\examples{
## A dog's preferences for six samples of food
mat <- matrix(c(0, 1, 1, 0, 1, 1,
                0, 0, 0, 1, 1, 0,
                0, 1, 0, 1, 1, 1,
                1, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 1,
                0, 1, 0, 1, 0, 0), 6, 6, byrow=TRUE)
circular(mat)  # moderate consistency
}
\keyword{models}
