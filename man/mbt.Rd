\name{mbt}
\alias{mbt}
\alias{print.mbt}
\title{Mallows-Bradley-Terry Model}
\description{
  Fits a Mallows-Bradley-Terry (MBT) model by maximum likelihood.
}
\usage{
mbt(data, bootstrap = FALSE, nsim = 1000, ...)
}
\arguments{
  \item{data}{a data frame, the first t columns containing the ranks, the
  (t + 1)th column containing the frequencies}
  \item{bootstrap}{logical. Return a parametric bootstrap p-value?}
  \item{nsim}{number of bootstrap replicates}
  \item{...}{further aguments passed to \code{simulate}}
}
\details{
  \code{mbt} provides a front end for \code{glm}.  See Critchlow and Fligner
  (1991) and Mallows (1957) for details.
}
\value{
  \item{coefficients}{a vector of parameter estimates (scale values)
    constrained to sum to unity}
  \item{goodness.of.fit}{the goodness of fit statistic including the
    likelihood ratio fitted vs. saturated model (-2logL), the degrees of
    freedom, the p-value of the corresponding chi-square distribution, and
    if \code{bootstrap} is \code{TRUE} the bootstrap p-value}
  \item{perm.idx}{the names of the non-zero frequency ranks}
  \item{y}{the vector of rank frequencies including zeros}
  \item{mbt.glm}{the output from a call to \code{glm}}
}
\author{Florian Wickelmaier}
\references{
  Critchlow, D.E., & Fligner, M.A. (1991).
  Paired comparison, triple comparison, and ranking experiments as generalized
  linear models, and their implementation in GLIM.
  \emph{Psychometrika}, \bold{56}, 517--533.
  \doi{10.1007/bf02294488}

  Mallows, C.L. (1957).
  Non-null ranking models. I.
  \emph{Biometrika}, \bold{44}, 114--130.
  \doi{10.1093/biomet/44.1-2.114}
}
\seealso{
  \code{\link{tartness}}, \code{\link{glm}}.
}
\examples{
data(tartness)        # tartness rankings of salad dressings (Vargo, 1989)
mbt(tartness, bootstrap=TRUE, nsim=500)  # fit Mallows-Bradley-Terry model
}
\keyword{models}
