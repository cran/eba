\name{boot}
\alias{boot}
\alias{eba.boot}
\title{Bootstrap for Elimination-By-Aspects Models}
\description{
  Performs a bootstrap by resampling the individual data matrices.
}
\usage{
  boot(D, R = 100, A = 1:I, s = rep(1/J, J))
}
\arguments{
  \item{D}{a 3d array consisting of the individual paired
    comparison matrices}
  \item{R}{the number of bootstrap samples}
  \item{A}{a list of vectors consisting of the stimulus aspects;
    the default is \code{1:I}, where \code{I} is the number of stimuli}
  \item{s}{the starting vector with default \code{1/J} for all parameters,
    where \code{J} is the number of parameters}
}
\details{
  The bootstrap function \code{eba.boot} is called automatically by
  \code{boot}.
}
\value{
  \item{p}{the matrix of bootstrap vectors}
  \item{stat}{the matrix of bootstrap statistics, including parameter
    means, standard errors, and confidence limits}}
}
\seealso{
  \code{\link{OptiPt}}.
}
\examples{
data(pork)  # pork tasting data, 10 individual paired comparison matrices
eba <- OptiPt(apply(pork, 1:2, sum))  # fit a BTL model
b <- boot(pork, 200)  # resample 200 times

plot(eba$estimate, b$stat[,'mean'], log='xy')
abline(0, 1, lty=2)
}
\keyword{models}
