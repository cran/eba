\name{residuals.eba}
\alias{residuals.eba}
\title{Residuals for EBA Models}
\description{
  Computes deviance and Pearson residuals for \code{eba} objects.
}
\usage{residuals.eba(object, type = c("deviance", "pearson"), ...)}
\arguments{
  \item{object}{an object of class \code{eba}, typically the result of a
    call to \code{\link{eba}}}
  \item{type}{the type of residuals which should be returned; the
          alternatives are: \code{"deviance"} (default) and
          \code{"pearson"}}
  \item{...}{further arguments passed to or from other methods.
    None are used in this method.}
}
\details{
  Residuals are computed from the upper triangle of the paired-comparison
  matrix.
  
  See \code{help(residuals.glm)} for details.
}
\value{
  A vector of residuals having as many elements as pairs of stimuli.
}
\seealso{\code{\link{eba}}, \code{\link{residuals.glm}},
  \code{\link{plot.eba}}.}
\examples{
data(celebrities)  # absolute choice frequencies
btl <- eba(celebrities)  # fit BTL model
sum( resid(btl)^2 )  # deviance
sum( resid(btl, "pearson")^2 )  # Pearson chi2
}
\keyword{models}
