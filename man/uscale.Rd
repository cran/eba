\name{uscale}
\alias{uscale}
\title{Utility Scale of an EBA Choice Model}
\description{
  Extract the (normalized) utility scale for an elimination-by-aspects (EBA)
  model.
}
\usage{uscale(object, norm = "sum", log = FALSE)}
\arguments{
  \item{object}{an object of class \code{eba}, typically the result of a
    call to \code{\link{eba}}}
  \item{norm}{either \code{sum}, so the scale values sum to unity (default);
    or a number from 1 to number of stimuli, so this scale value becomes one;
    or \code{NULL} (no normalization)}
  \item{log}{should the log of the utility scale values be returned? Defaults
    to FALSE.}
}
\details{
  Each utility scale value is defined as the sum of aspect values (EBA model
  parameters) that characterize a given stimulus.  First these sums are
  computed for all stimuli, then normalization (if any) is applied.  As
  each type of normalization corresponds to a multiplication by a positive
  real number, the ratio between scale values remains constant.
}
\value{
  The (normalized) utility scale of the stimuli.
}
\seealso{
  \code{\link{eba}}, \code{\link{cov.u}}, \code{\link{wald.test}}.
}
\examples{
data(drugrisk)
A <- list(c(1), c(2,7), c(3,7), c(4,7,8), c(5,7,8), c(6,7,8))
eba1 <- eba(drugrisk[, , group = "male30"], A)  # EBA model

uscale(eba1)                    # sum-to-unity normalization
uscale(eba1, norm=1)            # u(alcohol) := 1
uscale(eba1, norm=5)            # u(heroine) := 1
uscale(eba1, norm=NULL)         # no normalization
uscale(eba1, norm=1, log=TRUE)  # log utility scale, log u(alcohol) := 0
}
\keyword{models}
