\encoding{utf8}
\name{locmodes}
\alias{locmodes}
\alias{plot.locmod}
\alias{print.locmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Location of modes and antimodes }
\description{
   Given a certain number of modes, this function provides the estimation of the location of modes and antimodes and their density value.
}
\usage{
locmodes(data,mod0=1,lowsup=-Inf,uppsup=Inf,n=2^15,tol=10^(-5),display=F,\dots)

\method{plot}{locmod}(x,addplot=NULL,xlab=NULL,ylab=NULL,addLegend=NULL,posLegend=NULL,\dots)

\method{print}{locmod}(x,digits=getOption("digits"), \dots)
}
\arguments{
  \item{data}{ Sample in which the critical bandwidth is computed.}
  \item{mod0}{ Number of modes for which the critical bandwidth is calculated. Default \code{mod0=1}.}
  \item{lowsup}{ Lower limit for the random variable support in the computation of the  critical bandwidth. Default is \code{-Inf}.}
  \item{uppsup}{ Upper limit for the random variable support in the computation of the critical bandwidth. Default is \code{Inf}.}
  \item{n}{ The number of equally spaced points at which the density is to be estimated. When n > 512, it is rounded up to a power of 2 as in the \code{density} function. Default \code{n=2^15}.}
  \item{tol}{ Accuracy requested in the computation of the critical bandwidth. Default value \code{tol=10^(-5)}.}
  \item{display}{ Logical, if \code{TRUE} the kernel density estimation is plotted adding the estimated location of the modes and the antimodes. Default is \code{FALSE}.}
    \item{\dots}{ Arguments to be passed to subsequent methods, \code{\link[graphics]{plot.default}} for the \code{plot} method and \code{\link[base]{formatC}} for the \code{print} method.}
    \item{x}{ An object inheriting from class \code{"locmod"}.}
    \item{addplot}{ Logical, if \code{TRUE} the plot is added to the current one. Default is \code{FALSE}.}
    \item{xlab}{ A title for the x axis. See \code{title}.}
    \item{ylab}{ A title for the y axis. See \code{title}.}
  \item{addLegend}{ Logical, if \code{TRUE} the legend is added in the plot. Default is \code{TRUE}.}
  \item{posLegend}{ The vector of two elements of coordinates to be used to position the legend. It can be specified by keyword as in the function \code{legend}. Default is \code{"topright"}.}
    \item{digits}{ Number of significant digits to use, see \code{\link[base]{formatC}}.}
}

\details{

Given a certain number of modes, \code{mod0}, with \code{locmodes} the estimation of the location of modes and antimodes, their density value and the corresponding critical bandwidth is provided. To obtain these estimates, the kernel density estimation with gaussian kernel and the critical bandwidth for \code{mod0} modes is employed. If the compact support is unknown, the critical bandwidth of Silverman (1981) is computed and if such a support is provided, then the one proposed by Hall and York (2001) is calculated. Note that when the support is unknown the critical bandwidth may create artificial modes in the tails.

Since a dichotomy method is employed for computing the critical bandwidth, the parameter \code{tol} is used to determine a stopping time in such a way that the error committed in the computation of the critical bandwidth is less than \code{tol}.

If \code{display=TRUE}, then the kernel density estimation using the critical bandwidth for \code{mod0} modes is plotted. Additionally, the estimated location of modes (dashed lines), antimodes (point lines) and support (solid lines) can be also plotted. If \code{addLegend=TRUE}, a legend (in the position \code{posLegend}) with this information is included.

The NAs will be automatically removed.
}
\value{
  A list with class \code{"locmod"} containing the following components:
  \item{locations}{ Vector with the estimated locations of modes (odd positions of the vector) and antimodes (even positions). }
  \item{fvalue}{ Vector with estimated density values at modes (odd positions of the vector) and antimodes (even positions).}
  \item{cbw}{ A list with class \code{"estmod"} which contains the critical bandwidth of the sample for \code{mod0} modes, see \code{\link{bw.crit}}.}
}
\references{

Jose Ameijeiras--Alonso, Rosa M. Crujeiras, Alberto Rodríguez--Casal (2018). Mode testing, critical bandwidth and excess mass, to appear in \emph{Test}.

Hall, P. and York, M. (2001).  On the calibration of Silverman's test for multimodality, \emph{Statistica Sinica}, \bold{11}, 515--536.

Silverman, B. W. (1981). Using kernel density estimates to investigate multimodality, \emph{Journal of the Royal Statistical Society. Series B}, \bold{43}, 97--99.

}

\author{ Jose Ameijeiras--Alonso, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\examples{
# Testing for unimodality
set.seed(2016)
data=rnorm(50)
modetest(data)
#There is no evidence to reject the null hypothesis of unimodality
#Estimated location of the mode and its density value
locmodes(data)
\dontrun{
#Estimated locations of the five modes in the claw of Marron and Wand (1992)
library(nor1mix)
set.seed(2016)
n<-200
data<-nor1mix::rnorMix(n,MW.nm10)
#Adding the plot of the estimated locations
locmodes(data,5,display=T)
#Assuming that the compact support is [-1.5,1.5]
locmodes(data,5,-1.5,1.5,display=T)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mode estimation }% __ONLY ONE__ keyword per line

