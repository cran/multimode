\encoding{utf8}
\name{chondrite}
\alias{chondrite}
\alias{chondritegg}
\alias{chondritel}
\docType{data}
\title{ Percentage of silica in chondrite meteors }
\description{
This dataset contains the percentage silica (in \%) in 22 chondrite meteors.
}
\usage{
data(chondrite)
}
\format{
  \code{chondrite} and \code{chondritegg} include the data provided in Good and Gaskins  (1980). In \code{chondrite}, the typo in the 9th observation was corrected while \code{chondritegg} includes the original dataset. In \code{chondritel} is included the scaled data used by Leonard (1978).
}
\details{
This is a classic example for determining the number of modes.
}
\source{
Obtained from the Table 2 of Good and Gaskins (1980).
}
\references{

Good, I. J. and Gaskins, R. A. (1980). Density estimation and bump-hunting by the penalized likelihood method exemplified by scattering and meteorite data. \emph{Journal of the American Statistical Association}, \bold{75}, 42--56.

Leonard, T. (1978). Density estimation, stochastic processes and prior information. \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, \bold{40}, 113--146.

}
\examples{
data("chondrite")
# SiZer between the critical bandwidths for one and six modes
sizer(chondrite,cbw1=1,cbw2=6)

}
\keyword{datasets}
