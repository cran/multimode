\encoding{utf8}
\name{bw.crit}
\alias{bw.crit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Critical bandwidth }
\description{
   This function computes the critical bandwidth for a specified number of modes.
}
\usage{
bw.crit(data,mod0=1,lowsup=-Inf,uppsup=Inf,n=2^15,tol=10^(-5),full.result=F)
}
\arguments{
  \item{data}{ Sample for computing the critical bandwidth.}
  \item{mod0}{ Number of modes for which the critical bandwidth is calculated. Default \code{mod0=1}.}
  \item{lowsup}{ Lower limit for the random variable support in the computation of the  critical bandwidth. Default is \code{-Inf}.}
  \item{uppsup}{ Upper limit for the random variable support in the computation of the critical bandwidth. Default is \code{Inf}.}
  \item{n}{ The number of equally spaced points at which the density is estimated. When n > 512, it is rounded up to a power of 2 as in the \code{\link[stats]{density}} function. Default \code{n=2^15}.}
  \item{tol}{ Accuracy requested in the computation of the critical bandwidth. Default value of \code{tol} is \code{10^(-5)}.}
    \item{full.result}{ If this argument is TRUE then it returns the full result list, see below. Default \code{full.result=FALSE}.}
  }
\details{
With \code{bw.crit} the critical bandwidth for the number of modes specified in \code{mod0} is calculated, e.g., the smallest bandwidth such that the kernel density estimator has at most \code{mod0} modes. If the compact support is unknown, the critical bandwidth introduced by Silverman (1981) is computed and if it is provided that one of Hall and York (2001) is calculated.

Since a dichotomy method is employed for computing the critical bandwidth, the parameter \code{tol} is used to determine a stopping time in such a way that the error committed in the computation of the critical bandwidth is less than \code{tol}.

The NAs will be automatically removed.
}
\value{
  Depending on \code{full.result} either a number, the critical bandwidth of the sample for \code{mod0} modes, or an object of class \code{"estmod"} which is a \code{list} containing the following components:
  \item{nmodes}{ The specified hypothesized value of the number of modes.}
  \item{sample.size}{ The number of non-missing observations in the sample used for computing the number of modes.}
  \item{bw}{ Value of the critical bandwidth test statistic.}
  \item{lowsup}{ Lower limit of the support where the critical bandwidth is computed.}
  \item{ippsup}{ Upper limit of the support where the critical bandwidth is computed.}
  \item{fnx}{ The \code{n} coordinates of the points where the density is estimated for computing the critical bandwidth.}
  \item{fny}{ The estimated density values.}
}
\references{

Hall, P. and York, M. (2001).  On the calibration of Silverman's test for multimodality, \emph{Statistica Sinica}, \bold{11}, 515--536.

Silverman, B. W. (1981). Using kernel density estimates to investigate multimodality, \emph{Journal of the Royal Statistical Society. Series B}, \bold{43}, 97--99.

}

\author{ Jose Ameijeiras--Alonso, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\examples{
# Critical bandwidth of Silverman (1981) for one mode.
set.seed(2016)
data=rnorm(50)
bw.crit(data)
# Critical bandwidth of Hall and York for two modes in the interval (-1.5,1.5).
set.seed(2016)
data=rnorm(50)
bw.crit(data,mod0=2,lowsup=-1.5,uppsup=1.5)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mode testing }% __ONLY ONE__ keyword per line


