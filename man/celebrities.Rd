\name{celebrities}
\alias{celebrities}
\title{Choice among Celebrities}
\description{
  This data set provides the absolute choice frequencies of
  234 subjects choosing between pairs of nine celebrities.
  L. B. Johnson (LBJ), Harold Wilson (HW), and Charles De Gaulle
  (CDG) are politicians; Johnny Unitas (JU), Carl Yastrzemski
  (CY), and A. J. Foyt (AJF) are athletes; Brigitte Bardot
  (BB), Elizabeth Taylor (ET), and Sophia Loren (SL) are female
  movie stars.  Subjects were instructed to choose the person
  with whom they would rather spend an hour of discussion.
  Row stimuli are chosen over column stimuli.
}
\usage{data(celebrities)}
\format{A square data frame with a diagonal of zeros.}
\source{
  Rumelhart, D.L., & Greeno, J.G. (1971).
  Similarity between stimuli: An experimental test of the
  Luce and Restle choice models.
  \emph{Journal of Mathematical Psychology}, \bold{8}, 370--381.
}
\examples{
data(celebrities)
celebrities["LBJ", "HW"]  # 159 subjects choose LBJ over HW
}
\keyword{datasets}
