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
}
\value{
  A matrix having \code{(nstimuli-1)*nstimuli/2} rows and
  \code{nstimuli} columns.
}
\seealso{\code{\link{eba}}, \code{\link{glm}}.}
\examples{
data(celebrities)  # absolute choice frequencies
btl <- eba(celebrities)  # fit BTL model using eba
summary(btl)

y1 <- t(celebrities)[lower.tri(celebrities)]
y0 <- celebrities[lower.tri(celebrities)]

# Fit BTL model using glm
btl.glm <- glm(cbind(y1, y0) ~ pcX(9) - 1, binomial)
summary(btl.glm)  # BTL model does not fit these data!

# Fit Thurstone Case V model using glm
tm.glm <- glm(cbind(y1, y0) ~ pcX(9) - 1, binomial(probit))
summary(tm.glm)  # neither fits the TM model
}
\keyword{models}
