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
    absolute frequencies}
}
\details{
  With erroneous data, the number of violations is only a descriptive
  index of the goodness of fit of a model. When the sample size is small
  (or too large), however, the violations might be the only reliable
  index for the model selection.
}
\value{
  A table displaying the number of violations of the weak, moderate,
  and strong stochastic transitivity, the number of tests, the
  error ratio (violations/tests), and the mean and maximum deviation.
  \item{pcm}{the square matrix of relative choice frequencies}
}
\seealso{\code{\link{OptiPt}}.}
\examples{
data(rugr)  # absolute choice frequencies
strans(rugr)  # report of the transitivity violations
strans(rugr)$pcm  # relative frequencies
}
\keyword{models}
