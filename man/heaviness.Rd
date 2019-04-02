\name{heaviness}
\alias{heaviness}
\title{Weights Judging Data}
\description{
  Beaver and Gokhale (1975) presented fifty subjects with all 20 ordered pairs
  of bottles filled with lead shot and asked them to choose the bottle that
  felt heavier.  The mass of the bottles was 90, 95, 100, 105, and 110 grams,
  respectively.  Choice frequencies were aggregated across subjects for the
  two within-pair presentation orders.
}
\usage{data(heaviness)}
\format{
  A 3d array consisting of two square matrices:
  \describe{
    \item{\code{heaviness[, , order = 1]}}{holds the choices where the row
      stimulus was presented first for each pair (in the upper triangle, and
      vice versa in the lower triangle).}
    \item{\code{heaviness[, , order = 2]}}{holds the choices where the column
      stimulus was presented first for each pair (in the upper triangle, and
      vice versa in the lower triangle).}
  }
}
\source{
  Beaver, R.J., & Gokhale, D.V. (1975).
  A model to incorporate within-pair order effects in paired comparisons.
  \emph{Communications in Statistics}, \bold{4}, 923--939.
  \doi{10.1080/03610927308827302}
}
\seealso{\code{\link{eba.order}} for a model that includes a within-pair
  order effect.
}
\examples{
data(heaviness)
## 6 subjects chose 90g over 100g, when 90g was presented first.
heaviness["90g", "100g", order=1]

## 44 subjects chose 100g over 90g, when 90g was presented first.
heaviness["100g", "90g", order=1]

## 14 subjects chose 90g over 100g, when 90g was presented second.
heaviness["90g", "100g", order=2]

## 36 subjects chose 100g over 90g, when 90g was presented second.
heaviness["100g", "90g", order=2]

## Bradley-Terry-Luce (BTL) model for each within-pair order
btl1 <- eba(heaviness[, , 1])
btl2 <- eba(heaviness[, , 2])

xval <- seq(90, 110, 5)
plot(uscale(btl1) ~ xval, type="n", log="y",
     xlab="Mass of lead shot filled bottle (in g)",
     ylab="Perceived heaviness (BTL model)",
     main="Weights judging experiment")
mtext("(Beaver and Gokhale, 1975)", line=.5)

arrows(xval, uscale(btl1) - 1.96*sqrt(diag(cov.u(btl1))),
       xval, uscale(btl1) + 1.96*sqrt(diag(cov.u(btl1))), .05, 90, 3, "gray")
arrows(xval, uscale(btl2) - 1.96*sqrt(diag(cov.u(btl2))),
       xval, uscale(btl2) + 1.96*sqrt(diag(cov.u(btl2))), .05, 90, 3, "gray")

points(uscale(btl1) ~ xval, type="b", pch=16)
points(uscale(btl2) ~ xval, type="b", lty=2, pch=21, bg="white")

text(90.3, .01,  "Lower weight judged first", adj=0)
text(90.3, .08, "Higher weight judged first", adj=0)
}
\keyword{datasets}
