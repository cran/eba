\name{thurstone}
\alias{thurstone}
\alias{print.thurstone}
\title{Thurstone-Mosteller Model (Case V)}
\description{
  Fits a Thurstone-Mosteller model (Case V) by maximum likelihood.
}
\usage{
thurstone(M)
}
\arguments{
  \item{M}{a square matrix or a data frame consisting of absolute choice
    frequencies; row stimuli are chosen over column stimuli}
}
\details{
  \code{thurstone} provides a front end for \code{glm}.

  See Critchlow & Fligner (1991) for more details.
}
\value{
  \item{estimate}{a vector of parameter estimates (scale values), first
    element is set to zero}
  \item{goodness.of.fit}{the goodness of fit statistic including the
    likelihood ratio fitted vs. saturated model (-2logL), the degrees of
    freedom, and the p-value of the corresponding chi-square distribution}
  \item{tm.glm}{the output from a call to \code{glm}}
}
\author{Florian Wickelmaier}
\references{
  Critchlow, D.E., & Fligner, M.A. (1991). Paired comparison, triple
  comparison, and ranking experiments as generalized linear models, and their
  implementation in GLIM.
  \emph{Psychometrika}, \bold{56}, 517--533.
}
\seealso{
  \code{\link{eba}}, \code{\link{strans}}, \code{\link{pcX}},
  \code{\link{kendall.u}}, \code{\link{circular}}, \code{\link{glm}}.
}
\examples{
## Taste data (David, 1988, p. 116)
dat <- matrix(c( 0,  3,  2,  2,
                12,  0, 11,  3,
                13,  4,  0,  5,
                13, 12, 10,  0), 4, 4, byrow=TRUE)
thurstone(dat)  # Thurstone-Mosteller model fits OK
}
\keyword{models}
