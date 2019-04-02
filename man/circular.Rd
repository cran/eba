\name{circular}
\alias{circular}
\alias{print.circular}
\title{Circular Triads (Intransitive Cycles)}
\description{
  Number of circular triads and coefficient of consistency.
}
\usage{circular(mat, alternative = c("two.sided", "less", "greater"),
         exact = NULL, correct = TRUE, simulate.p.value = FALSE,
         nsim = 2000)}
\arguments{
  \item{mat}{a square matrix or a data frame consisting of (individual)
    binary choice data; row stimuli are chosen over column stimuli}
  \item{alternative}{a character string specifying the alternative hypothesis,
    must be one of \code{"two.sided"} (default), \code{"less"} or
    \code{"greater"}}
  \item{exact}{a logical indicating whether an exact p-value should be
    computed}
  \item{correct}{a logical indicating whether to apply continuity correction
    in the chi-square approximation for the p-value}
  \item{simulate.p.value}{a logical indicating whether to compute p-values by
    Monte Carlo simulation}
  \item{nsim}{an integer specifying the number of replicates used in the Monte
    Carlo test}
}
\details{
  Kendall's coefficient of consistency,
  \deqn{zeta = 1 - T/T_{max},}
  lies between one (perfect consistency) and zero,
  where \code{T} is the observed number of circular triads,
  and the maximum possible number of circular triads is
  \eqn{T_{max} = n(n^2 - 4)/24}, if \eqn{n} is even, and
  \eqn{T_{max} = n(n^2 - 1)/24} else, and \eqn{n} is the
  number of stimuli or objects being judged. For details see Kendall and
  Babington Smith (1940) and David (1988).

  Kendall (1962) discusses a test of the hypothesis that the number of
  circular triads \code{T} is different (smaller or greater) than expected
  when choosing randomly. For small \eqn{n}, an exact p-value is computed,
  based on the exact distributions listed in Alway (1962) and in Kendall
  (1962).  Otherwise, an approximate chi-square test is computed.  In this
  test, the sampling distribution is measured from lower to higher values of
  \code{T}, so that the probability that \code{T} will be exceeded is the
  complement of the probability for \code{chi2}.  The chi-square approximation
  may be incorrect if \eqn{n < 8} and is only available for \eqn{n > 4}.
}
\value{
  \item{T}{number of circular triads}
  \item{T.max}{maximum possible number of circular triads}
  \item{T.exp}{expected number of circular triads \eqn{E(T)} when choices are
    totally random}
  \item{zeta}{Kendall's coefficient of consistency}
  \item{chi2, df, correct}{the chi-square statistic and degrees of freedom for
    the approximate test, and whether continuity correction has been applied}
  \item{p.value}{the p-value for the test (see Details)}
  \item{simulate.p.value, nsim}{whether the p-value is based on simulations,
    number of simulation runs}
}
\references{
  Alway, G.G. (1962).
  The distribution of the number of circular triads in paired comparisons.
  \emph{Biometrika}, \bold{49}, 265--269.
  \doi{10.1093/biomet/49.1-2.265}

  David, H. (1988).
  \emph{The method of paired comparisons}. 
  London: Griffin.

  Kendall, M.G. (1962).
  \emph{Rank correlation methods}.
  London: Griffin.

  Kendall, M.G., & Babington Smith, B. (1940).
  On the method of paired comparisons. 
  \emph{Biometrika}, \bold{31}, 324--345.
  \doi{10.1093/biomet/31.3-4.324}
}
\seealso{
  \code{\link{eba}}, \code{\link{strans}}, \code{\link{kendall.u}}.
}
\examples{
# A dog's preferences for six samples of food
# (Kendall and Babington Smith, 1940, p. 326)
dog <- matrix(c(0, 1, 1, 0, 1, 1,
                0, 0, 0, 1, 1, 0,
                0, 1, 0, 1, 1, 1,
                1, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 1,
                0, 1, 0, 1, 0, 0), 6, 6, byrow=TRUE)
dimnames(dog) <- setNames(rep(list(c("meat", "biscuit", "chocolate",
                                     "apple", "pear", "cheese")), 2),
                          c(">", "<"))
circular(dog, alternative="less")  # moderate consistency
subset(strans(dog)$violdf, !wst)   # list circular triads
}
\keyword{models}
