\name{residuals.eba}
\alias{residuals.eba}
\title{Residuals for EBA Models}
\description{
  Computes deviance and Pearson residuals for \code{eba} objects.
}
\usage{residuals.eba(object, type = c("deviance", "pearson"), ...)}
\arguments{
  \item{object}{an object of class \code{eba}, typically the result of a
    call to \code{\link{OptiPt}}}
  \item{type}{the type of residuals which should be returned; the
          alternatives are: \code{"deviance"} (default) and
          \code{"pearson"}}
  \item{...}{further arguments passed to or from other methods}
}
\details{
  Residuals are computed from the upper triangle of the paired-comparison
  matrix.
  
  See \code{help(residuals.glm)} for details.
}
\value{
  A vector of residuals having as many elements as pairs of stimuli.
}
\seealso{\code{\link{OptiPt}}, \code{\link{residuals.glm}},
  \code{\link{plot.eba}}.}
\examples{
data(celebrities)  # absolute choice frequencies
btl <- OptiPt(celebrities)
sum( resid(btl)^2 )  # Deviance
sum( resid(btl, "pearson")^2 )  # Pearson chi2
}
\keyword{models}
