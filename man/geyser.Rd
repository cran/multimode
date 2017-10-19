\encoding{utf8}
\name{geyser}
\alias{geyser}
\alias{geyserab}
\alias{geyserh}
\alias{geyserw}
\docType{data}
\title{ Waiting time between geyser eruptions }
\description{
This dataset contains the interval times between the starts of the geyser eruptions on the Old Faithful Geyser in Yellowstone National Park, Wyoming, USA.
}
\usage{
data(geyser)
data(geyserab)
data(geyserh)
data(geyserw)
}
\format{
  \code{geyser}, \code{geyserh} and \code{geyserw} include data collected in October 1980. \code{geyser} and \code{geyserh} were obtained from Table 3 of Härdle (2012), in \code{geyser} some repeated data are removed. \code{geyserw} was obtained from the Supplementary material of Weisberg (2005). \code{geyserab} include data collected in August 1985 from Table 1 in Azzalini and Bowman (1990).
}
\details{
This is a classic example for estimating the density.
}
\source{
Obtained from the original tables in Azzalini and Bowman (1990), Härdle (2012) and Supplementary material of Weisberg (2005).
}
\references{

Azzalini, A. and Bowman, A. W. (1990). A look at some data on the Old Faithful geyser. \emph{Applied Statistics}, \bold{39}, 357-365.

Härdle, W. (1991). \emph{Smoothing techniques: with implementation in S}. New York: Springer-Verlag.

Weisberg, S. (2005). \emph{Applied Linear Regression}. New York: Wiley.

}
\examples{
data("geyser")
# Kernel density estimation with two modes
locmodes(geyser,mod0=2,display=TRUE)

}
\keyword{datasets}
