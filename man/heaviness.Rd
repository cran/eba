\name{heaviness}
\alias{heaviness}
\title{Weights Judging Data}
\description{
  Fifty subjects were presented with all 20 ordered pairs of
  bottles filled with lead shot and asked to choose the
  bottle that felt heavier.  The weight of the bottles was
  90, 95, 100, 105, and 110 grams, respectively.  Choice
  frequencies were aggregated across subjects for the two
  within-pair presentation orders.
}
\usage{data(heaviness)}
\format{
  A 3d array consisting of two square matrices:
  \describe{
    \item{\code{heaviness[,,order=1]}}{holds the choices where the row
      stimulus was presented first for each pair (in the upper triangle, and
      vice versa in the lower triangle).}
    \item{\code{heaviness[,,order=2]}}{holds the choices where the column
      stimulus was presented first for each pair (in the upper triangle, and
      vice versa in the lower triangle).}
  }
}
\source{
  Beaver, R.J., & Gokhale, D.V. (1975).
  A model to incorporate within-pair order effects in paired comparisons.
  \emph{Communications in Statistics}, \bold{4}, 923--939.
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
}
\keyword{datasets}
