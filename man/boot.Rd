\name{boot}
\alias{boot}
\alias{eba.boot}
\title{Bootstrap for Elimination-by-aspects models}
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
    the default is 1:I, where I is the number of stimuli}
  \item{s}{the starting vector with default 1/J for all parameters,
    where J is the number of parameters}
}
\details{
  The bootstrap function eba.boot is called automatically by boot.
}
\value{
  \item{p}{the matrix of bootstrap vectors}
  \item{stat}{the matrix of bootstrap statistics, including parameter
    means, standard errors, and confidence limits}}
}
\author{Florian Wickelmaier \email{wickelmaier@web.de}}
\seealso{
  \code{\link{OptiPt}}.
}
\examples{
data(pork)  # pork tasting data, 10 individual paired comparison matrices
eba = OptiPt(apply(pork,1:2,sum))  # fit a BTL model
b = boot(pork,200)  # resample 200 times

plot(eba$estimate,b$stat[,'mean'],log='xy')
abline(0,1,lty=3)
}
\keyword{models}
