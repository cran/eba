\name{mbt}
\alias{mbt}
\alias{print.mbt}
\title{Mallows-Bradley-Terry Model}
\description{
  Fits a Mallows-Bradley-Terry (MBT) model by maximum likelihood.
}
\usage{
mbt(data)
}
\arguments{
  \item{data}{a data frame, the first t columns containing the ranks, the
  (t + 1)th column containing the frequencies.}
}
\details{
  \code{mbt} provides a front end for \code{glm}.

  See Critchlow & Fligner (1991) for more details.
}
\value{
  \item{coefficients}{a vector of parameter estimates (scale values)
    constrained to sum to unity}
  \item{goodness.of.fit}{the goodness of fit statistic including the
    likelihood ratio fitted vs. saturated model (-2logL), the degrees of
    freedom, and the p-value of the corresponding chi-square distribution}
  \item{perm.idx}{a vector of indices of non-zero frequency ranks}
  \item{y}{the vector of rank frequencies including zeros}
  \item{mbt.glm}{the output from a call to \code{glm}}
}
\author{Florian Wickelmaier}
\references{
  Critchlow, D.E., & Fligner, M.A. (1991). Paired comparison, triple
  comparison, and ranking experiments as generalized linear models, and their
  implementation in GLIM.
  \emph{Psychometrika}, \bold{56}, 517--533.

  Mallows, C.L. (1957). Non-null ranking models. I. \emph{Biometrika},
  \bold{44}, 114--130.
}
\seealso{
  \code{\link{tartness}}, \code{\link{glm}}.
}
\examples{
data(tartness)  # tartness rankings of salad dressings (Vargo, 1989)
mbt(tartness)   # fit Mallows-Bradley-Terry model
}
\keyword{models}
