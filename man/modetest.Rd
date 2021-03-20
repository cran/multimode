\encoding{utf8}
\name{modetest}
\alias{modetest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Test for the number of modes }
\description{
   This function tests the number of modes.
}
\usage{
modetest(data,mod0=1,method="ACR",B=500,lowsup=-Inf,uppsup=Inf,
submethod=NULL,n=NULL,tol=NULL,tol2=NULL,gridsize=NULL,alpha=NULL,
nMC=NULL,BMC=NULL)
}
\arguments{
  \item{data}{ Sample to be tested.}
  \item{mod0}{ The maximum number of modes in the null hypothesis. Default \code{mod0=1} (unimodality vs. multimodality test).}
  \item{method}{ The method employed for testing the number of modes. Available methods are: SI (Silverman, 1981), HY (Hall and York, 2001), FM (Fisher and Marron, 2001), HH (Hartigan and Hartigan, 1985), CH (Cheng and Hall,1998), ACR (Ameijeiras-Alonso \emph{et al.}, 2019). Default \code{method="ACR"}.}
  \item{B}{ Number of replicates used in the test. Default \code{B=500}.}
  \item{lowsup}{ Lower limit for the random variable support in the computation of the  critical bandwidth. Default is \code{-Inf}.}
  \item{uppsup}{ Upper limit for the random variable support in the computation of the critical bandwidth. Default is \code{Inf}.}
    \item{submethod}{ Different approaches when using \code{method} SI, HY or ACR. Available submethods are: 1 and 2, see Details below for more information. Default \code{submethod=1}, except when \code{mod0}>1, \code{method} is ACR and the sample size is greater than 200.}
  \item{n}{ The number of equally spaced points at which the density is to be estimated. When n > 512, it is rounded up to a power of 2 as for \code{\link[stats]{density}} function. Default is \code{n=2^10} when \code{method} is SI, HY or FM and \code{n=2^15} when the \code{method} is CH or ACR. Argument not used for other methods.}
  \item{tol}{ Accuracy for computing the critical bandwidth. Default \code{tol=10^(-5)}.}
  \item{tol2}{ Accuracy for integration of the calibration function in the \code{method} ACR when the support is known. Default \code{tol2=10^(-5)}.}
  \item{gridsize}{ When the approximated version of the excess mass is employed in \code{method} ACR, number of endpoints at which the \eqn{C_m(\lambda)} sets are estimated (first element) and number of possible values of \eqn{\lambda} (second element). Default is \code{gridsize=c(20,20)}.}
  \item{alpha}{ Significance level employed for testing unimodality when method 1 of Hall and York (2001) is used. Default \code{alpha=0.05}.}
  \item{nMC}{ Number of Monte Carlo replicates used to approximate the p-value in the method 2 of Hall and York (2001). Default \code{nMC=100}.}
  \item{BMC}{ Number of bootstrap replicas used for computing the p-value in each Monte Carlo replicate of the Hall and York (2001) method 2. Default \code{BMC=100}.}
  }
\details{
The number of modes for the underlying density of a sample given by \code{data} can be tested with \code{modetest}. The null hypothesis states that the sample has \code{mod0} modes, and the alternative hypothesis is if it has more modes. The test used for calculating the p-value is specified in \code{method}. All the available proposals require bootstrap or Monte Carlo resamples (number specified in \code{B}).

Except when the support is employed, the typical usages are

\preformatted{
modetest(data,...)
modetest(data,mod0=1,method="ACR",B=500,...)
}

Since a dichotomy algorithm is employed for computing the critical bandwidth (methods SI, HY, FM, ACR), the parameter \code{tol} is used to determine a stopping time in such a way that the error committed in the computation of the critical bandwidth is less than \code{tol}.

The sample data can be perturbed in the methods using the excess mass or the dip statistic (HH, CH and ACR) in order to avoid important effects on the computation of the test statistic. In general, the exact excess mass/dip value is employed, but also its approximated version can be used by setting \code{submethod=2} in \code{method} ACR. See \code{\link{excessmass}} Details for more information.

Typical usages are

\preformatted{
modetest(data,mod0=1,method="ACR",B=500,submethod=1,n=NULL,tol=NULL)
modetest(data,mod0=1,method="ACR",B=500,submethod=2,n=NULL,tol=NULL,
gridsize=NULL)
}

When employing SI method, two ways of computing the resamples are available. If \code{submethod=1}, the resamples are generated from the rescaled bootstrap resamples as proposed by Silverman (1981). If \code{submethod=2}, as in \code{method} HY, the resamples are generated from the distribution that is associated to the kernel density estimation with the critical bandwidth.

Typical usage is

\preformatted{
modetest(data,mod0=1,method="SI",B=500,submethod=NULL,n=NULL,tol=NULL)
}

If a compact support containing a mode is known, it can be used to compute the Hall and York (2001) critical bandwidth. Note that in the case of the test proposed by Hall and York (2001), this support should be known, unless the support of the density function is bounded. For their proposal, two methods are implemented. \code{submethod} 1 is an asymptotic correction of Silverman (1981) test based on the limiting distribution of the test statistic. When \code{submethod=1}, the significance level must be previously determined with \code{alpha}. \code{submethod} 2 is based on Monte Carlo techniques. For this reason, when \code{submethod=2}, the number of replicates (\code{nMC}) and the number of bootstrap replicates (\code{BMC}) used for computing the p-value in each Monte Carlo replicate are needed.

Typical usages are

\preformatted{
modetest(data,method="HY",lowsup=-1.5,uppsup=1.5,...)
modetest(data,method="HY",B=500,lowsup=-1.5,uppsup=1.5,submethod=1,
n=NULL,tol=NULL,alpha=NULL)
modetest(data,method="HY",B=500,lowsup=-1.5,uppsup=1.5,submethod=2,
n=NULL,tol=NULL,nMC=NULL,BMC=NULL)
}

A modification of the proposal of Ameijeiras-Alonso \emph{et al.} (2019) can be also applied, by setting \code{method=ACR} and including a known compact support for detecting the modes. The parameter \code{tol2} is the accuracy required in the integration of the calibration function. For more information, see Ameijeiras-Alonso \emph{et al.} (2019), the default approach, when the support is unknown, is given in Section 2.3 and, when it is provided, the approach shown in the Appendix B is employed.

Typical usage is

\preformatted{
modetest(data,mod0=1,method="ACR",B=500,lowsup=-1.5,uppsup=1.5,
submethod=NULL,n=NULL,tol=NULL,tol2=NULL,gridsize=NULL)
}


The NAs will be automatically removed.
}


\value{
A list with class \code{"htest"} containing the following components:
  \item{p.value}{ P-value obtained after applying the test.}
  \item{statistic}{ Value of the test statistic. Critical bandwidth if the \code{method} is SI or HY; Cramer-von Mises statistic if the \code{method} is FM; the dip statistic if the \code{method} is HH, and the excess mass when the \code{method} is CH or ACR.}
  \item{null.value}{ The specified hypothesized value of the number of modes.}
  \item{alternative}{ A character string describing the alternative hypothesis, which is always \code{"greater"}.}
  \item{method}{ A character string indicating what type of multimodality test was performed.}
  \item{sample.size}{ The number of non-missing observations in the sample used for the hypothesis test.}
  \item{data.name}{ A character string giving the name of the data.}
  \item{bad.obs}{ The number of missing values that were removed from the data object prior to performing the hypothesis test.}
  }
\references{

Ameijeiras-Alonso, J., Crujeiras, R.M. and Rodríguez-Casal, A. (2019). Mode testing, critical bandwidth and excess mass, \emph{Test}, \bold{28}, 900--919.

Ameijeiras-Alonso, J., Crujeiras, R.M. and Rodríguez-Casal, A. (2021). multimode: An R Package for Mode Assessment, \emph{Journal of Statistical Software}, \bold{97}, 1--32.

Cheng, M. Y. and Hall, P. (1998). Calibrating the excess mass and dip tests of modality, \emph{Journal of the Royal Statistical Society. Series B}, \bold{60}, 579--589.

Fisher, N.I. and Marron, J. S. (2001). Mode testing via the excess mass estimate,
 \emph{Biometrika}, \bold{88}, 419--517.

Hall, P. and York, M. (2001).  On the calibration of Silverman's test for multimodality, \emph{Statistica Sinica}, \bold{11}, 515--536.

Hartigan, J. A. and Hartigan, P. M. (1985). The Dip Test of Unimodality,
 \emph{Journal of the American Statistical Association}, \bold{86}, 738--746.

Silverman, B. W. (1981). Using kernel density estimates to investigate multimodality, \emph{Journal of the Royal Statistical Society. Series B}, \bold{43}, 97--99.

}

\author{ Jose Ameijeiras-Alonso, Rosa M. Crujeiras and Alberto Rodríguez-Casal}

\examples{
# Testing for unimodality
data(geyser)
data=geyser
modetest(data)
\dontrun{
# Testing bimodality using B=100 bootstrap replicas
modetest(data,mod0=2,B=100)
#There is no evidence to reject the null hypothesis of bimodality
locmodes(data,mod0=2,display=TRUE)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mode testing }% __ONLY ONE__ keyword per line



