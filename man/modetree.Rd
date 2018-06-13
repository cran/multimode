\encoding{utf8}
\name{modetree}
\alias{modetree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mode tree }
\description{
   This function provides the mode tree.
}
\usage{
modetree(data,bws=NULL,gridsize=NULL,cbw1=NULL,cbw2=NULL,display=TRUE,
logbw=FALSE,\dots)
}
\arguments{
  \item{data}{ Sample in which the mode tree is computed.}
  \item{bws}{ Vector or range of bandwidths. If it is a vector of size two, then it is used a grid of bandwidths between the given values. Default lower bandwidth is twice the distance between the grid points used for estimating the density and upper bandwidth equal to the range of the data. Unless it is specified a vector of size greater than two, the number of bandwidths employed is equal to the second element of \code{gridsize}.}
  \item{gridsize}{ Number of equally spaced points at which the density is to be estimated (first element) and bandwidths used to compute the mode tree (second element). Default is \code{gridsize=c(512,151)}.}
  \item{cbw1}{ Number of modes for which the first critical bandwidth is calculated. This is the first bandwidth used to compute the mode tree when \code{bws} is not specified.}
  \item{cbw2}{ Number of modes for which the second critical bandwidth is calculated. This is the last bandwidth used to compute the mode tree when \code{bws} is not specified.}
  \item{display}{ Logical, if \code{TRUE} the mode tree plot is plotted. Default \code{TRUE}.}
  \item{logbw}{ Logical, if \code{TRUE} the plot displays and returns the log10 bandwidths. Default \code{logbw=FALSE}.}
  \item{\dots}{ Arguments to be passed to subsequent methods, \code{\link{plot.gtmod}} for the \code{plot} and \code{print} methods.}
}

\details{

The mode tree for the sample given in \code{data} is computed. For this calculation, a kernel density estimator with Gaussian kernel and bandwidths \code{bws} is used. When it is plotted, it shows with the continuous lines the estimated mode locations at each bandwidth. The horizontal dashed lines indicate the splitting of a mode in more modes.

The NAs will be automatically removed.
}
\value{
  A list with class \code{"gtmod"} containing the following components:
  \item{locations}{ Estimated mode locations for the bandwidths given in the row names.}
  \item{range.x}{ The data range, employed to represent the mode tree.}
  \item{range.bws}{ Employed bandwidths to compute the mode tree.}
  \item{logbw}{ Logical value indicating if the bandwidths are given in the log10 scale.}
  \item{sample.size}{ The number of non-missing observations in the sample used for obtaining the mode tree.}
  \item{call}{ The unevaluated expression, which consists of the named function applied to the given arguments.}
}
\references{

Minnotte, M. C. and Scott, D. W. (1993). The mode tree: A tool for visualization of nonparametric density features, \emph{Journal of Computational and Graphical Statistics}, \bold{2}, 51--68.

}


\author{ Jose Ameijeiras--Alonso, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\examples{
#Mode tree using a grid of bandwidths between 0.2 and 0.5
set.seed(2016)
data=rnorm(50)
modetree(data,bws=c(0.2,0.5))
\dontrun{
#Estimated locations of the five modes in the claw of Marron and Wand (1992)
library(nor1mix)
set.seed(2016)
n<-200
data<-nor1mix::rnorMix(n,MW.nm10)
#Mode tree between the critical bandwidths for 1 and 8 modes
modetree(data,cbw1=1,cbw2=8)
abline(v=1.5);abline(v=-1.5)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphical tools }% __ONLY ONE__ keyword per line

