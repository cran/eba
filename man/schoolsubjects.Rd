\name{schoolsubjects}
\alias{schoolsubjects}
\title{Preference for School Subjects}
\description{
  Two classes of children (ages 11 to 13) were asked to state their
  preferences with respect to certain school subjects. Each child was given a
  sheet on which were written the possible pairs of subjects and asked to
  underline the one preferred in each case (Kendall and Babington Smith,
  1940).
}
\usage{data(schoolsubjects)}
\format{
  A list containing two square matrices of aggregate choice frequencies
  (row entries are preferred over column entries):
  \describe{
    \item{\code{schoolsubjects[["boys"]]}}{holds the frequencies
      of 21 boys choosing among 13 school subjects:
      woodwork, gymnastics, art, science, history, geography, arithmetic,
      religion, English literature, commercial subjects, algebra,
      English grammar, geometry.}
    \item{\code{schoolsubjects[["girls"]]}}{holds the frequencies
      of 25 girls choosing among 11 school subjects:
      gymnastics, science, art, domestic science, history, arithmetic,
      geography, English literature, religion, algebra, English grammar.}
  }
}
\source{
  Kendall, M.G., & Babington Smith, B. (1940).
  On the method of paired comparisons.
  \emph{Biometrika}, \bold{31}, 324--345.
  \doi{10.1093/biomet/31.3-4.324}
}
\examples{
data(schoolsubjects)

m <- lapply(schoolsubjects, eba)  # Bradley-Terry-Luce (BTL) model

par(mfrow = 1:2)
dotchart(uscale(m$boys), main = "Boys' preferences",)
dotchart(uscale(m$girls), main = "Girls' preferences")
mtext("Utility scale value (BTL model)", outer = TRUE, side = 1,
      line = -2)
mtext("(Kendall and Babington Smith, 1940)", outer = TRUE, line = -4)
}
\keyword{datasets}

