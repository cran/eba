\name{pcX}
\alias{pcX}
\title{Paired-Comparison Design Matrix}
\description{
  Computes a paired-comparison design matrix.
}
\usage{pcX(nstimuli)}
\arguments{
  \item{nstimuli}{number of stimuli in the paired-comparison design}
}
\details{
  The design matrix can be used when fitting a Bradley-Terry-Luce (BTL)
  model or a Thurstone-Mosteller (TM) model by means of \code{\link{glm}}.

  See Critchlow & Fligner (1991) for more details.
}
\value{
  A matrix having \code{(nstimuli-1)*nstimuli/2} rows and
  \code{nstimuli} columns.
}
\references{
  Critchlow, D.E., & Fligner, M.A. (1991). Paired comparison, triple
  comparison, and ranking experiments as generalized linear models, and their
  implementation in GLIM.
  \emph{Psychometrika}, \bold{56}, 517--533.
}
\seealso{
  \code{\link{eba}}, \code{\link{thurstone}}, \code{\link{glm}},
  \code{\link{balanced.pcdesign}}.
}
\examples{
data(celebrities)        # absolute choice frequencies
btl <- eba(celebrities)  # fit BTL model using eba
summary(btl)

y1 <- t(celebrities)[lower.tri(celebrities)]
y0 <- celebrities[lower.tri(celebrities)]

## Fit BTL model using glm
btl.glm <- glm(cbind(y1, y0) ~ pcX(9) - 1, binomial)
summary(btl.glm)  # BTL model does not fit these data!

## Fit Thurstone Case V model using glm
tm.glm <- glm(cbind(y1, y0) ~ pcX(9) - 1, binomial(probit))
summary(tm.glm)   # neither fits the TM model
}
\keyword{models}
