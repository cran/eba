\name{weights}
\alias{weights}
\title{Weights Judging Data}
\description{
  Fifty subjects were presented with all 20 ordered pairs of
  bottles filled with lead shot and asked to choose the
  bottle that felt heavier.  The weight of the bottles was
  90, 95, 100, 105, and 110 grams, respectively.  Choice
  frequencies were aggregated across subjects for the two
  within-pair presentation orders.
}
\usage{data(weights)}
\format{
  A 3d array consisting of two square matrices:
  \describe{
    \item{\code{weights[,,"order1"]}}{holds the choices where the row
      stimulus was presented first for each pair (in the upper triangle, and
      vice versa in the lower triangle).}
    \item{\code{weights[,,"order2"]}}{holds the choices where the column
      stimulus was presented first for each pair (in the upper triangle, and
      vice versa in the lower triangle).}
  }
}
\source{
  Davidson, R.R., & Beaver, R.J. (1977).
  On extending the Bradley-Terry model to incorporate
  within-pair order effects.
  \emph{Biometrics}, \bold{33}, 693--702.
}
\examples{
data(weights)
# 6 subjects chose 90g over 100g, when 90g was presented first.
weights["90g", "100g", "order1"]

# 44 subjects chose 100g over 90g, when 90g was presented first.
weights["100g", "90g", "order1"]

# 14 subjects chose 90g over 100g, when 90g was presented second.
weights["90g", "100g", "order2"]

# 36 subjects chose 100g over 90g, when 90g was presented second.
weights["100g", "90g", "order2"]
}
\keyword{datasets}
