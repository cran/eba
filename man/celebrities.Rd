\name{celebrities}
\alias{celebrities}
\title{Choice among Celebrities}
\description{
  Rumelhart and Greeno (1971) presented 234 participants with pairs of names
  of nine celebrities: the politicians L. B. Johnson (LBJ), Harold Wilson
  (HW), and Charles De Gaulle (CDG); the athletes Johnny Unitas (JU), Carl
  Yastrzemski (CY), and A. J. Foyt (AJF); the actresses Brigitte Bardot (BB),
  Elizabeth Taylor (ET), and Sophia Loren (SL).  Participants were instructed
  to choose the person with whom they would rather spend an hour of
  discussion.
}
\usage{data(celebrities)}
\format{A square data frame containing the absolute choice frequencies and a
  diagonal of zeros; row stimuli are chosen over column stimuli.}
\source{
  Rumelhart, D.L., & Greeno, J.G. (1971).
  Similarity between stimuli: An experimental test of the
  Luce and Restle choice models.
  \emph{Journal of Mathematical Psychology}, \bold{8}, 370--381.
  \doi{10.1016/0022-2496(71)90038-1}
}
\examples{
data(celebrities)
celebrities["LBJ", "HW"]  # 159 participants chose Johnson over Wilson
}
\keyword{datasets}
