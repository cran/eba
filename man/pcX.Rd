\name{pcX}
\alias{pcX}
\title{Paired-Comparison Design Matrix}
\description{
  Computes a paired-comparison design matrix.
}
\usage{pcX(nstimuli, omitRef = TRUE)}
\arguments{
  \item{nstimuli}{number of stimuli in the paired-comparison design}
  \item{omitRef}{logical, if \code{TRUE} (default), the first column
    corresponding to the reference category is omitted}
}
\details{
  The design matrix can be used when fitting a Bradley-Terry-Luce (BTL)
  model or a Thurstone-Mosteller (TM) model by means of \code{\link{glm}}
  or \code{\link{lm}}.

  See Critchlow and Fligner (1991) for more details.
}
\value{
  A matrix having \code{(nstimuli - 1)*nstimuli/2} rows and
  \code{nstimuli - 1} columns (if the reference category is omitted).
}
\references{
  Critchlow, D.E., & Fligner, M.A. (1991).
  Paired comparison, triple comparison, and ranking experiments as generalized
  linear models, and their implementation in GLIM.
  \emph{Psychometrika}, \bold{56}, 517--533.
  \doi{10.1007/bf02294488}
}
\seealso{
  \code{\link{eba}}, \code{\link{thurstone}}, \code{\link{glm}},
  \code{\link{balanced.pcdesign}}, \code{\link{linear2btl}}.
}
\examples{
data(drugrisk)               # absolute choice frequencies
btl <- eba(drugrisk[, , 1])  # fit Bradley-Terry-Luce model using eba
summary(btl)

y1 <- t(drugrisk[, , 1])[lower.tri(drugrisk[, , 1])]
y0 <-   drugrisk[, , 1][ lower.tri(drugrisk[, , 1])]

## Fit Bradley-Terry-Luce model using glm
btl.glm <- glm(cbind(y1, y0) ~ 0 + pcX(6), binomial)
summary(btl.glm)

## Fit Thurstone Case V model using glm
tm.glm <- glm(cbind(y1, y0) ~ 0 + pcX(6), binomial(probit))
summary(tm.glm)
}
\keyword{models}
