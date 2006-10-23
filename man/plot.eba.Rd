\name{plot.eba}
\alias{plot.eba}
\title{Diagnostic Plot for EBA Models}
\description{
  Plots EBA residuals against fitted values.
}
\usage{plot.eba(x, ...)}
\arguments{
  \item{x}{an object of class \code{eba}, typically the result of a
    call to \code{\link{eba}}}
  \item{...}{further arguments passed to or from other methods.
    None are used in this method.}
}
\details{
  The deviance residuals are plotted against the predicted choice
  probabilities for the upper triangle of the paired-comparison
  matrix.
}
\seealso{\code{\link{eba}}, \code{\link{residuals.eba}}.}
\examples{
data(celebrities)  # absolute choice frequencies
btl <- eba(celebrities)  # fit BTL model
plot(btl)  # residuals versus fitted
}
\keyword{models}
