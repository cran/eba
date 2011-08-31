\name{residuals.eba}
\alias{residuals.eba}
\title{Residuals for EBA Models}
\description{
  Computes deviance and Pearson residuals for \code{eba} objects.
}
\usage{
\method{residuals}{eba}(object, type = c("deviance", "pearson"), \dots)
}
\arguments{
  \item{object}{an object of class \code{eba}, typically the result of a
    call to \code{\link{eba}}}
  \item{type}{the type of residuals which should be returned; the
          alternatives are: \code{"deviance"} (default) and
          \code{"pearson"}}
  \item{\dots}{further arguments passed to or from other methods.
    None are used in this method.}
}
\details{
  Residuals are computed from the upper triangle of the paired-comparison
  matrix.
  
  See \code{\link{residuals.glm}} for details.
}
\value{
  A vector of residuals having as many elements as pairs of stimuli.
}
\seealso{\code{\link{eba}}, \code{\link{residuals.glm}},
  \code{\link{plot.eba}}.}
\examples{
data(celebrities)                # absolute choice frequencies
btl1 <- eba(celebrities)         # fit Bradley-Terry-Luce model
sum( resid(btl1)^2 )             # deviance G2
sum( resid(btl1, "pearson")^2 )  # Pearson X2
}
\keyword{models}
