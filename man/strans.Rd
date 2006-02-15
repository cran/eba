\name{strans}
\alias{strans}
\alias{print.strans}
\title{Stochastic Transitivity}
\description{
  Checks the weak, moderate, and strong stochastic transitivity.
}
\usage{strans(M)}
\arguments{
  \item{M}{a square matrix or a data frame consisting of
    absolute choice frequencies}
}
\details{
  The weak (WST), moderate (MST), and strong (SST) stochastic
  transitivity hold for a set of choice probabilities \eqn{P},
  whenever if \eqn{P_{ij} \ge 0.5} and \eqn{P_{jk} \ge 0.5}, then

  \eqn{P_{ik} \ge 0.5} (WST),

  \eqn{P_{ik} \ge \min(P_{ij}, P_{jk})} (MST),

  \eqn{P_{ik} \ge \max(P_{ij}, P_{jk})} (SST).
}
\value{
  A table displaying the number of violations of the weak, moderate,
  and strong stochastic transitivity, the number of tests, the
  error ratio (violations/tests), and the mean and maximum deviation.
  \item{weak}{number of violations of WST}
  \item{moderate}{number of violations of MST}
  \item{strong}{number of violations of SST}
  \item{n.tests}{number of transitivity tests performed}
  \item{wst.violations}{a vector of WST violations}
  \item{mst.violations}{a vector of MST violations}
  \item{sst.violations}{a vector of SST violations}
  \item{pcm}{the square matrix of relative choice frequencies}
}
\seealso{\code{\link{OptiPt}}.}
\examples{
data(celebrities)  # absolute choice frequencies
strans(celebrities)  # WST and MST hold, but not SST
strans(celebrities)$pcm  # relative frequencies
}
\keyword{models}
