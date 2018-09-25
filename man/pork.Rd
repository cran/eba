\name{pork}
\alias{pork}
\title{Pork Tasting Data}
\description{
  Bradley and Terry (1952) provide the individual choice matrices of two
  judges choosing between pairs of three samples of pork meet.  The pigs had
  been fed on either corn (C), corn plus peanut supplement (Cp), or corn plus
  a large peanut supplement (CP).  Each judge does five repetitions.  The data
  are stored in a 3d array, the first five matrices of which correspond to the
  five repetitions of the first judge, the last five to the repetitions of the
  second judge.  Row stimuli are chosen (preferred) over column stimuli.
}
\usage{data(pork)}
\format{A 3d array consisting of ten square matrices.}
\source{
  Bradley, R.A., & Terry, M.E. (1952).
  Rank analysis of incomplete block designs.
  I. The method of paired comparisons.
  \emph{Biometrika}, \bold{39}, 324--345.
  \doi{10.1093/biomet/39.3-4.324}
}
\examples{
data(pork)
apply(pork, 1:2, sum)  # aggregate choice frequencies
}
\keyword{datasets}
