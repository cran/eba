\name{winetaste}
\alias{winetaste}
\alias{ambilight}
\alias{redwines}

\title{Wine Tasting Data}

\description{
  Paired comparison judgments for two wine tasting studies:
  \code{ambilight} includes the results of a study on the effect of ambient
  lighting on the flavor of wine; \code{redwines} includes judgments on the
  sensory quality of red wines.
}

\usage{data("winetaste")}

\format{
\code{ambilight} A data frame containing 230 observations on 10 variables:
\describe{
  \item{preference, fruitiness, spiciness, sweetness}{Paired comparison of
    class \code{\link[psychotools]{paircomp}}; judgments for one of the
    6 ordered pairs of the blue, red, and white lighting conditions.}
  \item{age}{subject age}
  \item{gender}{factor, subject gender}
  \item{sensesmell}{self-rating of sense of smell and taste.}
  \item{likewine}{self-rating of general liking of wine.}
  \item{drinkwine}{factor, frequency of wine consumption.}
  \item{redwhite}{factor, preference for red or white wine.}
}
\code{redwines} A data frame containing 11 observations on 7 variables:
\describe{
  \item{bitterness, fruitiness, sourness, roundness, preference}{Paired
    comparison of class \code{\link[psychotools]{paircomp}}; judgments for
    all 10 pairs from 5 red wines:  Primitivo di Manduria, Cotes du Rhone,
    Bourgogne, Shiraz cuvee, Cabernet Sauvignon.}
  \item{best}{factor, Which of the five wines did you like best?}
  \item{worst}{factor, Which of the five wines did you like worst?}
}
}

\details{
  The \code{ambilight} data are from Experiment 3 in Oberfeld et al. (2009).
  The \code{redwines} data were collected among the members of the Sound
  Quality Research Unit (SQRU), Department of Acoustics, Aalborg University,
  Denmark, in 2004. Details of the red wines are available as an attribute of
  the \code{preference} variable (see Examples).
}

\references{
  Oberfeld, D., Hecht, H., Allendorf, U., & Wickelmaier, F. (2009).
  Ambient lighting modifies the flavor of wine.
  \emph{Journal of Sensory Studies},
  \bold{24}(6), 797--832.
  \doi{10.1111/j.1745-459X.2009.00239.x}
}

\seealso{
  \code{\link{eba}}, \code{\link{eba.order}},
  \code{\link[psychotools]{paircomp}}.
}

\examples{
requireNamespace("psychotools")
data(winetaste)

## No effect of ambient lighting on flavor (Oberfeld et al., 2009)

m <- lapply(ambilight[, c("preference", "fruitiness",
                          "spiciness", "sweetness")],
            function(x) eba.order(summary(x, pcmatrix = TRUE)))
lapply(m, summary)

u <- sapply(m, uscale, norm = 3)
dotchart(
  u, xlim = c(0.5, 2), pch = 16, panel.first = abline(v = 1, col = "gray"),
  main = "Ambient lighting and the flavor of wine",
  xlab = "Utility scale value (Davidson-Beaver model)"
)
ci <- sapply(m, function(x) 1.96 * sqrt(diag(cov.u(x))))
arrows(
  u - ci, c(16:18, 11:13, 6:8, 1:3), u + ci, c(16:18, 11:13, 6:8, 1:3),
.05, 90, 3)
mtext("Oberfeld et al. (2009)", line = 0.5)

## Sensory quality of red wines

psychotools::covariates(redwines$preference)  # details of the wines

m <- lapply(redwines[, c("bitterness", "fruitiness", "sourness",
                         "roundness", "preference")],
            function(x) eba(summary(x, pcmatrix = TRUE)))
lapply(m, summary)

u <- sapply(m, uscale)
dotchart(
  u[order(u[, "preference"]), ], log = "x",
  panel.first = abline(v = 1/5, col = "gray"),
  main = "SQRU red wine tasting",
  xlab = "Utility scale value (BTL model), choice proportion (+)"
)
points(as.vector(
  prop.table(table(redwines$best))[order(u[, "preference"])]
), 1:5, pch = 3)
}

\keyword{datasets}
