\name{OptiPt}
\alias{OptiPt}
\alias{summary.eba}
\alias{print.summary.eba}
\alias{print.eba}
\alias{L}
\title{Elimination-by-aspects models}
\description{
  Fits a probabilistic choice model by maximum likelihood estimation.
}
\usage{
  OptiPt(M, A = 1:I, s = rep(1/J, J))

  summary.eba(object, ...)
}
\arguments{
  \item{M}{a square matrix or a data frame consisting of
    absolute frequencies}
  \item{A}{a list of vectors consisting of the stimulus aspects;
    the default is 1:I, where I is the number of stimuli}
  \item{s}{the starting vector with default 1/J for all parameters,
    where J is the number of parameters}
  \item{object}{an object of class \code{eba}}
  \item{...}{additional arguments affecting the summary produced}
}
\details{
  A is usually a list of vectors, the first element of each being a
  number from 1 to I; additional elements specify the aspects shared
  by several stimuli.
  When fitting a BTL model, A reduces to 1:I, i.e. there is one aspect
  per stimulus.

  The maximum likelihood estimation of the parameter values is
  carried out by \code{nlm}. The Hessian matrix, however, is
  approximated by \code{fdHess{nlme}}.
}
\value{
  \item{estimate}{a vector of parameter estimates}
  \item{se}{a vector of the standard errors of the parameter estimates}
  \item{ci95}{the 95\%-confidence intervals for the parameter estimates}
  \item{logL.eba}{the log-likelihood of the fitted model}
  \item{logL.sat}{the log-likelihood of the saturated (binomial) model}
  \item{goodness.of.fit}{the goodness of fit statistic including the
    likelihood ratio (-2logL), the degrees of freedom, and the p-value
    of the corresponig chi2 distribution}
  \item{u.scale}{the u-scale of the stimuli; one scale value is defined
    as the sum of aspect values (parameters) that characterize a given
    stimulus}
  \item{hessian}{the Hessian matrix of the likelihood function}
  \item{cov.p}{the covariance matrix of the model parameters}
  \item{chi.alt}{the Pearson chi2 goodness of fit statistic}
  \item{fitted}{the fitted paired-comparison matrix}
}
\references{
  Wickelmaier, F., & Schmid, C. (2004). A Matlab function to estimate
  choice model parameters from paired-comparison data.
  \emph{Behavior Research Methods, Instruments, and Computers},
  \bold{36}, 29-40.

  Bradley, R.A. (1984). Paired comparisons: some basic procedures and
  examples. In P.R. Krishnaiah & P.K. Sen (eds.),
  \emph{Handbook of Statistics, Volume 4}. Amsterdam: Elsevier.

  Tversky, A. (1972). Elimination by aspects: A theory of choice.
  \emph{Psychological Review}, \bold{79}, 281-299.

  Tversky, A., & Sattath, S. (1979). Preference trees.
  \emph{Psychological Review}, \bold{86}, 542-573.
}
\author{Florian Wickelmaier \email{wickelmaier@web.de}}
\seealso{
  \code{\link{strans}}, \code{\link{cov.u}}, \code{\link{nlm}}.
}
\examples{
data(rugr)  # absolute choice frequencies
A = list(c(1,10),c(2,10),c(3,10),c(4,11),c(5,11),c(6,11),
         c(7,12),c(8,12),c(9,12))  # the structure of aspects
eba = OptiPt(rugr,A)  # Fit a preference tree

summary(eba)  # goodness of fit
ci = 1.96*sqrt(diag(cov.u(eba)))  # 95\%-ci of the preference scale

a = barplot(eba$u,ylim=c(0,.35))  # plot the scale
arrows(a,eba$u-ci,a,eba$u+ci,.05,90,3)  # error bars
}
\keyword{models}
