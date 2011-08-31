\name{plot.eba}
\alias{plot.eba}
\title{Diagnostic Plot for EBA Models}
\description{
  Plots elimination-by-aspects (EBA) model residuals against fitted values.
}
\usage{
\method{plot}{eba}(x, xlab = "Predicted choice probabilities",
  ylab = "Deviance residuals", \dots)
}
\arguments{
  \item{x}{an object of class \code{eba}, typically the result of a
    call to \code{\link{eba}}}
  \item{xlab, ylab, \dots}{graphical parameters passed to plot.}
}
\details{
  The deviance residuals are plotted against the predicted choice
  probabilities for the upper triangle of the paired-comparison
  matrix.
}
\seealso{\code{\link{eba}}, \code{\link{residuals.eba}}.}
\examples{
## Compare two choice models

data(celebrities)                     # absolute choice frequencies
btl1 <- eba(celebrities)              # fit Bradley-Terry-Luce model
A <- list(c(1,10), c(2,10), c(3,10),
          c(4,11), c(5,11), c(6,11),
          c(7,12), c(8,12), c(9,12))  # the structure of aspects
eba1 <- eba(celebrities, A)           # fit elimination-by-aspects model
anova(btl1, eba1)                     # model comparison based on likelihoods

par(mfrow = c(1,2))                   # residuals versus fitted values
plot(btl1, main = "BTL", ylim = c(-4, 4.5))  # BTL doesn't fit well
plot(eba1, main = "EBA", ylim = c(-4, 4.5))  # EBA fits better
}
\keyword{models}
