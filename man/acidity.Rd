\encoding{utf8}
\name{acidity}
\alias{acidity}
\docType{data}
\title{ Acid--neutralizing capacity }
\description{
This dataset, analyzed by Crawford (1994), contains the Acid--Neutralizing Capacity (ANC) measured in a sample of 155
lakes in North--Central Wisconsin (USA). ANC describes the capability of a lake to absorb acid, where low ANC values may lead to a loss of biological resources.
}
\usage{
data(acidity)
}
\format{
  \code{acidity} includes an acidity index of lakes in north--central Wisconsin on the log scale, in particular, it is provided the log(ANC+50) as in Crawford (1994).
}
\details{
This is a classic example for determining the number of modes.
}
\source{
Obtained from the Supplementary material of Richardson and Green (1997), available in \url{http://www.stats.bris.ac.uk/~peter/mixdata}.
}
\references{

Crawford (1994). An application of the Laplace method to finite mixture distributions. \emph{Journal of the American Statistical Association}, \bold{89}, 259--267.

Richardson, S. and Green, P. J. (1997). On Bayesian analysis of mixtures with unknown number of components (with discussion). \emph{Journal of the Royal Statistical Society, Series B}, \bold{59}, 731--792.

}
\examples{
data("acidity")
# Kernel density estimation with two modes and SiZer
locmodes(acidity,mod0=2,display=TRUE,xlab="log(ANC+50)")
sizer(acidity,bws=c(0.1,1),xlab="log(ANC+50)")
}
\keyword{datasets}
