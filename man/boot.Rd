\name{boot}
\alias{boot}
\alias{eba.boot}
\alias{eba.boot.constrained}
\title{Bootstrap for Elimination-by-Aspects (EBA) Models}
\description{
  Performs a bootstrap by resampling the individual data matrices.
}
\usage{
boot(D, R = 100, A = 1:I, s = rep(1/J, J), constrained = TRUE)
}
\arguments{
  \item{D}{either a 3d array consisting of the individual paired
    comparison matrices or an object of class
    \code{\link[psychotools]{paircomp}}}
  \item{R}{the number of bootstrap samples}
  \item{A}{a list of vectors consisting of the stimulus aspects;
    the default is \code{1:I}, where \code{I} is the number of stimuli}
  \item{s}{the starting vector with default \code{1/J} for all parameters,
    where \code{J} is the number of parameters}
  \item{constrained}{logical, if TRUE (default), parameters are constrained
    to be positive}
}
\details{
  The bootstrap functions \code{eba.boot.constrained} and \code{eba.boot}
  are automatically called by \code{boot}.

  The code is experimental and may change in the future.
}
\value{
  \item{p}{the matrix of bootstrap vectors}
  \item{stat}{the matrix of bootstrap statistics, including parameter
    means, standard errors, and confidence limits}
}
\seealso{
  \code{\link{eba}}, \code{\link{simulate.eba}},
  \code{\link[psychotools]{paircomp}}.
}
\examples{
data(pork)  # pork tasting data, 10 individual paired comparison matrices
m <- eba(apply(pork, 1:2, sum))  # fit Bradley-Terry-Luce model
b <- boot(pork, R = 200)         # resample 200 times

plot(coef(m), b$stat[, "mean"], log = "xy")
abline(0, 1, lty = 2)
}
\keyword{models}
