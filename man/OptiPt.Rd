\name{eba}
\alias{OptiPt}
\alias{eba}
\alias{summary.eba}
\alias{anova.eba}
\alias{vcov.eba}
\alias{print.summary.eba}
\alias{print.eba}
\alias{L}
\alias{L.constrained}
\title{Elimination-By-Aspects (EBA) Models}
\description{
  Fits a (multi-attribute) probabilistic choice model by maximum likelihood.
}
\usage{
eba(M, A = 1:I, s = rep(1/J, J), constrained = TRUE)

OptiPt(M, A = 1:I, s = rep(1/J, J), constrained = TRUE)

\method{summary}{eba}(object, \dots)

\method{anova}{eba}(object, \dots, test = c("Chisq", "none"))
}
\arguments{
  \item{M}{a square matrix or a data frame consisting of absolute choice
    frequencies; row stimuli are chosen over column stimuli}
  \item{A}{a list of vectors consisting of the stimulus aspects;
    the default is \code{1:I}, where \code{I} is the number of stimuli}
  \item{s}{the starting vector with default \code{1/J} for all parameters,
    where \code{J} is the number of parameters}
  \item{constrained}{logical, if TRUE (default), parameters are constrained
    to be positive}
  \item{object}{an object of class \code{eba}, typically the result of a
    call to \code{\link{eba}}}
  \item{test}{should the p-values of the chi-square distributions be
    reported?}
  \item{\dots}{additional arguments; none are used in the summary method;
    in the anova method they refer to additional objects of class \code{eba}.}
}
\details{
  \code{eba} is a wrapper function for \code{OptiPt}.  Both functions
  can be used interchangeably.

  The probabilistic choice models that can be fitted to paired-comparison
  data are the Bradley-Terry-Luce (BTL) model, preference tree (Pretree)
  models, and elimination-by-aspects (EBA) models, the former being special
  cases of the latter.

  \code{A} is usually a list of vectors, the first element of each being
  a number from 1 to \code{I}; additional elements specify the aspects shared
  by several stimuli.  \code{A} must have as many elements as there are
  stimuli.  When fitting a BTL model, \code{A} reduces to \code{1:I} (the
  default), i.e. there is only one aspect per stimulus.

  The maximum likelihood estimation of the parameter values is
  carried out by \code{nlm}.  The Hessian matrix, however, is
  approximated by \code{fdHess{nlme}}.  The likelihood functions
  \code{L.constrained} and \code{L} are called automatically.

  See \code{\link{group.test}} for details on the likelihood ratio
  tests reported by \code{summary.eba}.
}
\value{
  \item{estimate}{a vector of parameter estimates}
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
  \item{y1}{the data vector of the upper triangle matrix}
  \item{y0}{the data vector of the lower triangle matrix}
  \item{n}{the number of observations per pair (\code{y1 + y0})}
  \item{mu}{the predicted choice probabilities for the upper triangle}
}
\author{Florian Wickelmaier}
\references{
  Wickelmaier, F., & Schmid, C. (2004). A Matlab function to estimate
  choice model parameters from paired-comparison data.
  \emph{Behavior Research Methods, Instruments, and Computers},
  \bold{36}, 29--40.

  Bradley, R.A. (1984). Paired comparisons: some basic procedures and
  examples. In P.R. Krishnaiah & P.K. Sen (eds.),
  \emph{Handbook of Statistics, Volume 4}. Amsterdam: Elsevier.

  Tversky, A. (1972). Elimination by aspects: a theory of choice.
  \emph{Psychological Review}, \bold{79}, 281--299.

  Tversky, A., & Sattath, S. (1979). Preference trees.
  \emph{Psychological Review}, \bold{86}, 542--573.
}
\seealso{
  \code{\link{strans}}, \code{\link{cov.u}}, \code{\link{wald.test}},
  \code{\link{group.test}}, \code{\link{plot.eba}},
  \code{\link{residuals.eba}}, \code{\link{logLik.eba}},
  \code{\link{kendall.u}}, \code{\link{circular}}, \code{\link{thurstone}},
  \code{\link{nlm}}.
}
\examples{
data(celebrities)                     # absolute choice frequencies
btl1 <- eba(celebrities)              # fit a Bradley-Terry-Luce model
A <- list(c(1,10), c(2,10), c(3,10), c(4,11), c(5,11), c(6,11),
          c(7,12), c(8,12), c(9,12))  # the structure of aspects
eba1 <- eba(celebrities, A)           # fit a preference tree model

summary(eba1)      # goodness of fit
plot(eba1)         # residuals versus predicted values
anova(btl1, eba1)  # model comparison based on likelihoods

ci <- 1.96 * sqrt(diag(cov.u(eba1)))    # 95\% CI for preference scale values
a <- barplot(eba1$u, ylim = c(0, .35))  # plot the scale
arrows(a, eba1$u-ci, a, eba1$u+ci, .05, 90, 3)  # error bars
}
\keyword{models}
