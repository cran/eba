\name{trineq}
\alias{trineq}
\alias{print.trineq}
\title{Trinary Inequality}
\description{
  Checks if binary choice probabilities fulfill the trinary inequality.
}
\usage{trineq(M, A = 1:I)}
\arguments{
  \item{M}{a square matrix or a data frame consisting of absolute choice
    frequencies; row stimuli are chosen over column stimuli}
  \item{A}{a list of vectors consisting of the stimulus aspects;
    the default is \code{1:I}, where \code{I} is the number of stimuli}
}
\details{
  For any triple of stimuli \eqn{x, y, z}, the trinary inequality states
  that, if \eqn{P(x, y) > 1/2} and \eqn{(xy)z}, then
  \deqn{R(x, y, z) > 1,}
  where \eqn{R(x, y, z) = R(x, y) R(y, z) R(z, x)},
  \eqn{R(x, y) = P(x, y)/P(y, x)}, and \eqn{(xy)z} denotes that \eqn{x} and
  \eqn{y} share at least one aspect that \eqn{z} does not have
  (Tversky & Sattath, 1979, p. 554).

  \code{\link{inclusion.rule}} checks if a family of aspect sets is
  representable by a tree.
}
\value{
  Results checking the trinary inequality.
  \item{n}{number of tests of the trinary inequality}
  \item{prop}{proportion of triples confirming the trinary inequality}
  \item{quant}{quantiles of \eqn{R(x, y, z)}}
  \item{n.tests}{number of transitivity tests performed}
  \item{chkdf}{data frame reporting \eqn{R(x, y, z)} for each triple where
    \eqn{P(x, y) > 1/2} and \eqn{(xy)z}}
}
\references{
  Tversky, A., & Sattath, S. (1979). Preference trees.
  \emph{Psychological Review}, \bold{86}, 542--573.
}
\seealso{\code{\link{eba}}, \code{\link{inclusion.rule}},
  \code{\link{strans}}.}
\examples{
data(celebrities)             # absolute choice frequencies
A <- list(c(1,10), c(2,10), c(3,10),
          c(4,11), c(5,11), c(6,11),
          c(7,12), c(8,12), c(9,12))  # the structure of aspects
trineq(celebrities, A)        # check trinary inequality for tree A
trineq(celebrities, A)$chkdf  # trinary inequality for each triple
}
\keyword{models}
