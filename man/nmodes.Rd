\encoding{utf8}
\name{nmodes}
\alias{nmodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Number of modes }
\description{
   This function computes the number of modes in a kernel density estimator using the Gaussian kernel and a given bandwidth parameter.
}
\usage{
nmodes(data,bw,lowsup=-Inf,uppsup=Inf,n=2^15,full.result=F)
}
\arguments{
  \item{data}{ Sample for computing a kernel density estimator and determine the number of modes.}
  \item{bw}{ Bandwidth parameter for kernel density estimation.}
  \item{lowsup}{ Lower limit for the random variable support. Just the number of modes greater than \code{lowsup} are taken into account. Default is \code{-Inf}.}
  \item{uppsup}{ Upper limit for the random variable support. Just the number of modes greater than \code{lowsup} are taken into account. Default is \code{Inf}.}
  \item{n}{ The number of equally spaced points at which the density is to be estimated. When n > 512, it is rounded up to a power of 2 as in the \code{\link[stats]{density}} function. Default \code{n=2^15}.}
    \item{full.result}{ If this argument is TRUE then it returns the full result list, see below. Default \code{full.result=FALSE}.}
  }
\details{
The number of modes in the interval provided by \code{lowsup} and \code{uppsup} is computed. For this calculation, a kernel density estimator with Gaussian kernel and bandwidth \code{bw} is used.

The NAs will be automatically removed.
}
\value{
  Depending on \code{full.result} either a number, the number of modes for the bandwidth provided in \code{bw}, or an object of class \code{"estmod"} which is a \code{list} containing the following components:
  \item{nmodes}{ The number of modes for the bandwidth provided in \code{bw}.}
  \item{sample.size}{ The number of non-missing observations in the sample used for computing the number of modes.}
  \item{bw}{ Employed bandwidth for  kernel density estimation.}
  \item{lowsup}{ Lower limit of the support where the number of modes are computed.}
  \item{ippsup}{ Upper limit of the support where the number of modes are computed.}
  \item{fnx}{ The \code{n} coordinates of the points where the density is estimated for computing the number of modes.}
  \item{fny}{ The estimated density values.}
}

\author{ Jose Ameijeiras--Alonso, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\examples{
# Number of modes in the interval (-1.5,1.5), using the bandwidth 0.5.
set.seed(2016)
data=rnorm(50)
nmodes(data,0.5,-1.5,1.5)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mode estimation }% __ONLY ONE__ keyword per line

