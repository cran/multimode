\encoding{utf8}
\name{enzyme}
\alias{enzyme}
\docType{data}
\title{ Blood enzymatic activity }
\description{
This dataset concerns the distribution of enzymatic activity in the blood, for an enzyme involved in the metabolism of carcinogenic substances.
}
\usage{
data(enzyme)
}
\format{
  \code{enzyme} includes the values of the urinary metabolic ratio of 5--acetylamino--6--formylamino--3--methyluracil to 1--methylxanthine (AFMU/1X) after oral administration of caffeine.
}
\details{
This is a classic example for determining the number of modes.
}
\source{
Obtained from the Supplementary material of Richardson and Green (1997), available in \url{http://www.stats.bris.ac.uk/~peter/mixdata}.
}
\references{

Bechtel, Y. C., Bonaiti-Pellie, C., Poisson, N., Magnette, J. and Bechtel, P. R. (1993). A population and family study N--acetyltransferase using caffeine urinary metabolites. \emph{Clinical Pharmacology and Therapeutics}, \bold{54}, 134--141.

Richardson, S. and Green, P. J. (1997). On Bayesian analysis of mixtures with unknown number of components (with discussion). \emph{Journal of the Royal Statistical Society, Series B}, \bold{59}, 731--792.

}
\examples{
data("enzyme")
# It can be seen that there are two groups in this dataset
# Via exploratory tools
sizer(enzyme,bws=c(0.03,1))
modetree(enzyme,bws=c(0.02,1),logbw=TRUE,addplot=TRUE,col.lines="white")
#Via mode testing
modetest(enzyme)
\dontrun{
modetest(enzyme,mod0=2)
}
}
\keyword{datasets}
