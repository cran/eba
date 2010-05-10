\name{group.test}
\alias{group.test}
\alias{print.group.test}
\title{Group Effects in EBA Models}
\description{
  Tests for group effects in EBA models.
}
\usage{group.test(groups, A = 1:I, s = rep(1/J, J), constrained = TRUE)}
\arguments{
  \item{groups}{a 3d array containing one aggregate choice matrix
    per group}
  \item{A}{a list of vectors consisting of the stimulus aspects;
    the default is \code{1:I}, where \code{I} is the number of stimuli}
  \item{s}{the starting vector with default \code{1/J} for all parameters,
    where \code{J} is the number of parameters}
  \item{constrained}{logical, if TRUE (default), EBA parameters are
    constrained to be positive}
}
\details{
  The five tests are all based on likelihood ratios. 

  \code{Overall} compares a 1-parameter Poisson model to a saturated
  Poisson model, thereby testing the equality of the frequencies in each cell
  of the array.  This test corresponds to simultaneously testing for a null
  effect of (1) the context induced by a given pair, (2) the grouping factor,
  (3) the stimuli, and (4) the imbalance between pairs.  The deviances of
  the remaining tests sum to the total deviance associated with the overall
  test.

  \code{EBA.g} tests an EBA group model against a saturated
  binomial group model, which corresponds to a goodness of fit test of
  the EBA group model.

  \code{Group} tests an EBA model having its parameters restricted to
  be equal across groups (single set of
  parameters) against the EBA group model allowing its parameters to vary
  freely across groups (one set of parameters per
  group); this corresponds to testing for group differences.

  \code{Effect} tests an indifference model (where all choice probabilities
  equal 0.5) against the restricted EBA model,
  which corresponds to testing for a stimulus effect.

  \code{Imbalance} tests for differences in the number of observations
  per pair by comparing the average sample size (1-parameter Poisson model)
  to the actual sample sizes (saturated Poisson model).
  
  See Duineveld, Arents, & King (2000) for further details, and Choisel &
  Wickelmaier (2007) for an application.
}
\value{
  \item{tests}{a table displaying the likelihood ratio test statistics}
}
\references{
  Duineveld, C.A.A., Arents, P., & King, B.M. (2000). Log-linear modelling of
  paired comparison data from consumer tests.
  \emph{Food Quality and Preference}, \bold{11}, 63--70.

  Choisel, S., & Wickelmaier, F. (2007). Evaluation of multichannel
  reproduced sound: Scaling auditory attributes underlying listener
  preference.
  \emph{Journal of the Acoustical Society of America}, \bold{121},
  388--400.
}
\seealso{\code{\link{eba}}, \code{\link{wald.test}}.}
\examples{
data(pork)          # Is there a difference between Judge 1 and Judge 2?
groups <- array(c(apply(pork[,,1:5], 1:2, sum),
                  apply(pork[,,6:10], 1:2, sum)), c(3,3,2))
group.test(groups)  # Yes, there is.

data(drugrisk)   # Do younger and older males judge risk of drugs differently?
A2 <- list (c(1), c(2,7), c(3,7), c(4,7,8), c(5,7,8), c(6,7,8))
group.test(drugrisk[,,3:4], A2)  # Yes.
}
\keyword{models}
