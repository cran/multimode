\encoding{utf8}
\name{galaxy}
\alias{galaxy}
\alias{galaxyp}
\alias{galaxyrg}
\docType{data}
\title{ Velocities of galaxies diverging away from our own galaxy }
\description{
This dataset contains the velocities in km/sec of different galaxies from the unfilled survey of the Corona Borealis region.
}
\usage{
data(galaxy)
}
\format{
  \code{galaxy} includes the original measures in Roeder (1990). \code{galaxyrg} includes the data provided in the Supplementary material of Richardson and Green (1997), where the velocities are divided by 1000 and the 78th observation was replaced by 26690 km/sec. \code{galaxyp} add a measurement of 5607 km/sec included in Postman \emph{et al.} (1986).
}
\details{
This is a classic example for determining the number of modes.
}
\source{
Obtained from the Table 1 of Postman \emph{et al.} (1986), Table 1 of Roeder (1990) and the Supplementary material of Richardson and Green (1997), available in \url{http://www.stats.bris.ac.uk/~peter/mixdata}.
}
\references{

Postman, M., Huchra, J. P. and Geller, M. J. (1986). Probes of large-scale structures in the Corona Borealis region. \emph{Astronomical Journal}, \bold{92}, 1238--1247.

Richardson, S. and Green, P. J. (1997). On Bayesian analysis of mixtures with unknown number of components (with discussion). \emph{Journal of the Royal Statistical Society, Series B}, \bold{59}, 731--792.

Roeder, K. (1990). Density estimation with confidence sets exemplified by superclusters and voids in the galaxies. \emph{Journal of the American Statistical Association}, \bold{85}, 617--624.

}
\examples{
data("galaxy")
# Mode tree between the critical bandwidths for one and six modes
modetree(galaxy,cbw1=1,cbw2=6)

}
\keyword{datasets}
