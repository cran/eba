\name{linear2btl}
\alias{linear2btl}
\title{Linear Coefficients to Bradley-Terry-Luce (BTL) Estimates}
\description{
  Transforms linear model coefficients to Bradley-Terry-Luce (BTL) model
  parameter estimates.
}
\usage{linear2btl(object, order = FALSE)}
\arguments{
  \item{object}{an object of class \code{glm} or \code{lm} specifying a BTL
    model}
  \item{order}{logical, does the model include an order effect? Defaults to
    FALSE}
}
\details{
  The design matrix used by \code{glm} or \code{lm} usually results from
  a call to \code{\link{pcX}}. It is assumed that the reference category is
  the first level.

  The covariance matrix is estimated by employing the delta method.

  See Imrey, Johnson, & Koch (1976) for more details.
}
\value{
  \item{btl.parameters}{a matrix; the first column holds the BTL parameter
    estimates, the second column the approximate standard errors}
  \item{cova}{the approximate covariance matrix of the BTL parameter
    estimates}
  \item{linear.coefs}{a vector of the original linear coefficients as returned
    by \code{glm} or \code{lm}}
}
\references{
  Imrey, P.B., Johnson, W.D., & Koch, G.G. (1976). An incomplete contingency
  table approach to paired-comparison experiments.
  \emph{Journal of the American Statistical Association}, \bold{71}, 614--623.
}
\seealso{
  \code{\link{eba}}, \code{\link{eba.order}}, \code{\link{glm}},
  \code{\link{pcX}}.
}
\examples{
data(drugrisk)
y1 <- t(drugrisk[,,1])[lower.tri(drugrisk[,,1])]
y0 <-   drugrisk[,,1][ lower.tri(drugrisk[,,1])]

## Fit BTL model using glm (maximum likelihood)
btl.glm <- glm(cbind(y1, y0) ~ 0 + pcX(6), binomial)
linear2btl(btl.glm)

## Fit BTL model using lm (weighted least squares)
btl.lm <- lm(log(y1/y0) ~ 0 + pcX(6), weights=y1*y0/(y1 + y0))
linear2btl(btl.lm)
}
\keyword{models}
