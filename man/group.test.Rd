\name{group.test}
\alias{group.test}
\alias{print.group.test}
\title{Group Effects in EBA Models}
\description{
  Tests for group effects in EBA models.
}
\usage{group.test(groups, A = 1:I, s = rep(1/J, J))}
\arguments{
  \item{groups}{a 3d array including one aggregated choice matrix
    per group}
  \item{A}{a list of vectors consisting of the stimulus aspects;
    the default is \code{1:I}, where \code{I} is the number of stimuli}
  \item{s}{the starting vector with default \code{1/J} for all parameters,
    where \code{J} is the number of parameters}
}
\details{
  The five tests are all based on likelihood ratios. 
  \code{Overall} tests a 1-parameter poisson model against a saturated
  poisson model.
  \code{EBA.g} tests an EBA group model against a saturated
  binomial group model, which corresponds to a goodness of fit test of
  the EBA group model.
  \code{Group} tests an EBA model for the pooled data against the EBA
  group model, which corresponds to testing for group differences.
  \code{Effect} tests an indifference model against the pooled EBA model.
  \code{Imbalance} tests for differences in the number of observations
  per pair by comparing the average sample size (1-parameter poisson model)
  to the actual sample sizes (saturated poisson model).
  
  See Duineveld, Arents & King (2000) for details.
}
\value{
  \item{tests}{a table displaying the likelihood ratio test statistics}
}
\references{
  Duineveld, C.A.A., Arents, P., & King, B.M. (2000). Log-linear modelling of
  paired comparison data from consumer tests.
  \emph{Food Quality and Preference}, \bold{11}, 63-70.
}
\seealso{\code{\link{OptiPt}}, \code{\link{wald.test}}.}
\examples{
data(pork)  # Is there a difference between Judge 1 and Judge 2?
groups <- array(c(apply(pork[,,1:5], 1:2, sum),
                  apply(pork[,,6:10], 1:2, sum)), c(3,3,2))
group.test(groups)
}
\keyword{models}
