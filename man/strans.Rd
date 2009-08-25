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

  See Suppes, Krantz, Luce, & Tversky (1989/2007, chap. 17) for an
  introduction to the representation of choice probabilities.

  If WST holds, a permutation of the indices of the matrix exists such that
  the proportions in the upper triangular matrix are \eqn{\ge 0.5}. This
  re-arranged matrix is stored in \code{pcm}. If WST does not hold, cells in
  the upper triangular matrix that are smaller than 0.5 are replaced by 0.5.
  The deviance resulting from this restriction is reported in \code{wst.fit}.

  The approximate likelihood ratio test for significance of the WST violations
  is according to Tversky (1969); a more exact test of WST is suggested by
  Iverson & Falmagne (1985).
}
\value{
  A table displaying the number of violations of the weak, moderate,
  and strong stochastic transitivity, the number of tests, the
  error ratio (violations/tests), and the mean and maximum deviation
  from the minimum probability for which the corresponding transitivity
  would hold.
  \item{weak}{number of violations of WST}
  \item{moderate}{number of violations of MST}
  \item{strong}{number of violations of SST}
  \item{n.tests}{number of transitivity tests performed}
  \item{wst.violations}{a vector containing
    \eqn{0.5 - P_{ik}} for all triples that violate WST}
  \item{mst.violations}{a vector containing
    \eqn{\min(P_{ij}, P_{jk}) - P_{ik}} for all triples that violate MST}
  \item{sst.violations}{a vector containing
    \eqn{\max(P_{ij}, P_{jk}) - P_{ik}} for all triples that violate SST}
  \item{pcm}{the permuted square matrix of relative choice frequencies}
  \item{ranking}{the ranking of the objects, which corresponds to the colnames
    of pcm}
  \item{chkdf}{data.frame reporting the choice proportions for each triple in
    each permutation}
  \item{violdf}{data.frame reporting for each triple which type of
    transitivity holds or does not hold}
  \item{wst.fit}{likelihood ratio test of WST (see details)}
  \item{wst.mat}{restricted matrix that satisfies WST}
}
\references{
  Iverson, G., & Falmagne, J.-C. (1985). Statistical issues in measurement.
  \emph{Mathematical Social Sciences}, \bold{10}, 131--153.

  Suppes, P., Krantz, D.H., Luce, R.D., & Tversky, A. (1989/2007).
  \emph{Foundations of measurement. Volume II}. Mineola, N.Y.: Dover
  Publications.

  Tversky, A. (1969). Intransitivity of preferences.
  \emph{Psychological Review}, \bold{76}, 31--48.
}
\seealso{\code{\link{eba}}, \code{\link{circular}}, \code{\link{kendall.u}}.}
\examples{
data(celebrities)           # absolute choice frequencies
strans(celebrities)         # WST and MST hold, but not SST
strans(celebrities)$pcm     # re-ordered relative frequencies
strans(celebrities)$violdf  # log transitivity violations
}
\keyword{models}
