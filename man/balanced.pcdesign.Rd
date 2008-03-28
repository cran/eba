\name{balanced.pcdesign}
\alias{balanced.pcdesign}
\title{Balanced Paired-Comparison Design}
\description{
  Creates a (completely) balanced paired-comparison design.
}
\usage{balanced.pcdesign(nstimuli)}
\arguments{
  \item{nstimuli}{number of stimuli in the paired-comparison design}
}
\details{
  When \code{nstimuli} is odd, the presentation order is completely balanced,
  that is any given stimulus appears an equal number of times as the first
  and second member of a pair. When \code{nstimuli} is even, the presentation
  order is balanced as much as possible.

  In order to balance the within-pair presentation order over a sample of
  subjects, subjects should be equally assigned to \code{listA} and
  \code{listB}.

  Pairs should be re-randomized for each subject.
}
\value{
  \item{pairs}{a character array holding the balanced pairs; see David (1988)
    for details how it is constructed}
  \item{listA}{the vector pairs in the original within-pair order}
  \item{listB}{the vector of pairs in the inverted within-pair order}
}
\references{
  David, H. (1988). \emph{The method of paired comparisons}. 
  London: Griffin.
}
\seealso{
  \code{\link{pcX}}, \code{\link{eba}}.
}
\examples{
## Create balanced design for 6 stimuli
bp <- balanced.pcdesign(6)

## Replicate each with-pair order 10 times and re-randomize
cbind(replicate(10, sample(bp$listA)), replicate(10, sample(bp$listB)))
}
\keyword{models}
