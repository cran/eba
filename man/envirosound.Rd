\name{envirosound}
\alias{envirosound}

\title{Auditory Unpleasantness of Environmental Sound}

\description{
  Zimmer et al. (2004) investigated the auditory unpleasantness of
  twelve short binaural recordings (Johannsen and Prante, 2001); recordings
  were presented via headphones to 74 participants.
}

\usage{data(envirosound)}

\format{
  A data frame containing 74 observations on 2 variables:
  \describe{
  \item{unpleasantness}{paired comparison of class \code{\link{paircomp}};
    judgments for all 66 paired comparisons from 12 recordings:
    circular saw, stadium, dentist's drill, waterfall, ship's horn, stone in
    well, typewriter, hooves, fan, howling wind, tyre on gravel, wasp.}
    \item{rt}{median response time.}
  }
}

\details{
  Details of the recordings, including psychoacoustic metrics, are available
  as an attribute of the \code{unpleasantness} variable (see Examples).
}

\source{
  Zimmer, K., Ellermeier, W., & Schmid, C. (2004).
  Using probabilistic choice models to investigate auditory unpleasantness.
  \emph{Acta Acustica united with Acustica}, \bold{90}(6), 1019--1028.
}

\references{
  Johannsen, K., & Prante, H.U. (2001).
  Environmental sounds for psychoacoustic testing.
  \emph{Acta Acustica united with Acustica}, \bold{87}(2), 290--293.
}

\seealso{\code{\link{eba}}, \code{\link{strans}}, \code{\link{paircomp}}.}

\examples{
requireNamespace("psychotools")
data(envirosound)
set.seed(1019)

## Choice-model representation of unpleasantness
mat <- summary(envirosound$unpleasantness, pcmatrix = TRUE)
strans(mat)
btl1 <- eba(mat)
eba1 <- eba(mat, A = list(c(1, 13), c(2, 13), c(3, 13), c(4, 13),
                          c(5, 13), c(6, 13), c(7, 13), c(8, 13),
                          c(9, 13), c(10, 13), c(11, 13), 12))
eba2 <- eba(mat, A = list(c(1, 13), c(2, 13), c(3, 13), c(4, 13),
                          c(5, 13), c(6, 13), c(7, 13, 14), c(8, 13, 14),
                          c(9, 13, 14), c(10, 13, 14), c(11, 13, 14), 12),
            s = runif(14))
anova(btl1, eba1, eba2)

sounds <- psychotools::covariates(envirosound$unpleasantness)
sounds$u <- 10 * uscale(eba2, norm = 9)  # u(fan) := 10

plot(magnitude ~ u, sounds, log = "x", type = "n",
     xlab = "Indirect scaling (EBA model)",
     ylab = "Direct magnitude estimation",
     main = "Auditory unpleasantness of environmental sound")
mtext("(Zimmer et al., 2004)", line = 0.5)
abline(lm(magnitude ~ log10(u), sounds))
text(magnitude ~ u, sounds, labels = abbreviate(rownames(sounds), 4))

## Predicting unpleasantness from psychoacoustic metrics
summary(
  lm(log(u) ~ scale(sharpness, scale = FALSE) +
              scale(roughness, scale = FALSE):I(loudness.5 > 27),
     sounds[-12, ])  # w/o wasp
)
}
\keyword{datasets}
