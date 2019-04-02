\name{inclusion.rule}
\alias{inclusion.rule}
\title{Inclusion Rule}
\description{
  Checks if a family of sets fulfills the inclusion rule.
}
\usage{inclusion.rule(A)}
\arguments{
  \item{A}{a list of vectors consisting of the stimulus aspects of an
    elimination-by-aspects model}
}
\details{
  The inclusion rule is necessary and sufficient for a tree structure on the
  aspect sets:

  Structure theorem. A family \eqn{\{x' | x \in T\}} of aspect sets is
  representable by a tree iff either \eqn{x' \cap y' \supset x' \cap z'} or
  \eqn{x' \cap z' \supset x' \cap y'} for all \eqn{x, y, z} in \eqn{T}.
  (Tversky and Sattath, 1979, p. 546)
}
\value{
  Either \code{TRUE} if the inclusion rule holds for \code{A}, or \code{FALSE}
  otherwise.
}
\references{
  Tversky, A., & Sattath, S. (1979).
  Preference trees.
  \emph{Psychological Review}, \bold{86}, 542--573.
  \doi{10.1037/0033-295X.86.6.542}
}
\seealso{\code{\link{eba}}, \code{\link{trineq}}, \code{\link{strans}}.}
\examples{
A <- list(c(1, 5), c(2, 5), c(3, 6), c(4, 6))  # tree
inclusion.rule(A)

B <- list(c(1, 5), c(2, 5, 6), c(3, 6), c(4, 6))  # lattice
inclusion.rule(B)
}
\keyword{models}
