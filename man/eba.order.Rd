\name{eba.order}
\alias{eba.order}
\alias{summary.eba.order}
\alias{print.eba.order}
\alias{print.summary.eba.order}
\alias{L.order}
\alias{L.constrained.order}
\alias{logLik.eba.order}
\alias{plot.eba.order}
\alias{residuals.eba.order}
\title{Elimination-By-Aspects (EBA) Models with order-effect}
\description{
  Fits a (multi-attribute) probabilistic choice model that
  accounts for the effect of the presentation order within a pair.
}
\usage{
eba.order(M1, M2, A = 1:I, s = c(rep(1/J, J), 1), constrained = TRUE)

\method{summary}{eba.order}(object, \dots)
}
\arguments{
  \item{M1, M2}{two square matrices or data frames consisting of absolute
    choice frequencies in both within-pair orders; row stimuli are chosen
    over column stimuli}
  \item{A}{see \code{\link{eba}}}
  \item{s}{the starting vector with default \code{1/J} for all \code{J}
    aspect parameters, and \code{1} for the order effect}
  \item{constrained}{see \code{\link{eba}}}
  \item{object}{an object of class \code{eba.order}, typically the result
    of a call to \code{\link{eba.order}}}
  \item{\dots}{additional arguments}
}
\details{
  The choice models include a single multiplicative order effect,
  \code{order}, that is constant for all pairs (see Davidson & Beaver,
  1977). An order effect < 1 (> 1) indicates a bias in favor of the
  first (second) interval.

  See \code{\link{eba}} for choice models without order effect.

  Several likelihood ratio tests are performed (see also
  \code{\link{summary.eba}}).

  \code{EBA.order} tests an order-effect EBA model against a saturated
  binomial model; this corresponds to a goodness of fit test of
  the former model.

  \code{Order} tests an EBA model with an order effect constrained
  to 1 against an unconstrained order-effect EBA model; this corresponds
  to a test of the order effect.

  \code{Effect} tests an order-effect indifference model (where all
  scale values are equal, but the order effect is free) against
  the order-effect EBA model; this corresponds to testing for a
  stimulus effect; \code{order0} is the estimate of the former model.

  Wickelmaier & Choisel (2006) describe a model that generalizes the
  Davidson-Beaver model and allows for an order effect in Pretree and
  EBA models.
}
\value{
  \item{estimate}{a vector of parameter estimates, the last component holds
    the order-effect estimate}
  \item{se}{a vector of standard errors of the parameter estimates}
  \item{ci95}{a vector of 95\% confidence intervals for the parameters}
  \item{logL.eba}{the log-likelihood of the fitted model}
  \item{logL.sat}{the log-likelihood of the saturated (binomial) model}
  \item{goodness.of.fit}{the goodness of fit statistic including the
    likelihood ratio fitted vs. saturated model (-2logL), the degrees of
    freedom, and the p-value of the corresponding chi-square distribution}
  \item{u.scale}{the u-scale of the stimuli; one u-scale value is defined
    as the sum of aspect values (parameters) that characterize a given
    stimulus}
  \item{hessian}{the Hessian matrix of the likelihood function}
  \item{cov.p}{the covariance matrix of the model parameters}
  \item{chi.alt}{the Pearson chi-square goodness of fit statistic}
  \item{fitted}{the fitted paired-comparison matrix}
  \item{y1}{the data vector of the upper triangle matrices}
  \item{y0}{the data vector of the lower triangle matrices}
  \item{n}{the number of observations per pair (\code{y1 + y0})}
  \item{mu}{the predicted choice probabilities for the upper triangles}
  \item{M1, M2}{the data matrices}
}
\author{Florian Wickelmaier}
\references{
  Davidson, R.R., & Beaver, R.J. (1977).
  On extending the Bradley-Terry model to incorporate
  within-pair order effects.
  \emph{Biometrics}, \bold{33}, 693--702.

  Wickelmaier, F., & Choisel, S. (2006).
  Modeling within-pair order effects in paired-comparison judgments.
  In D.E. Kornbrot, R.M. Msetfi, & A.W. MacRae (eds.),
  \emph{Fechner Day 2006. Proceedings of the 22nd Annual Meeting of
  the International Society for Psychophysics} (p. 89--94).
  St. Albans, UK: The ISP.
}
\seealso{
  \code{\link{eba}}, \code{\link{group.test}}, \code{\link{plot.eba}},
  \code{\link{residuals.eba}}, \code{\link{logLik.eba}}.
}
\examples{
data(heaviness)  # weights judging data
ebao1 <- eba.order(heaviness[,,1], heaviness[,,2])  # Davidson-Beaver model
summary(ebao1)   # goodness of fit
plot(ebao1)      # residuals versus predicted values
}
\keyword{models}