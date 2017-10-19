\encoding{utf8}
\name{modetest}
\alias{modetest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Test for the number of modes }
\description{
   This function tests the number of modes.
}
\usage{
modetest(data,mod0=1,method="NP",B=500,full.result=FALSE,lowsup=-Inf,uppsup=Inf,
n=NULL,tol=NULL,tol2=NULL,methodhy=NULL,alpha=NULL,nMC=NULL,BMC=NULL)
}
\arguments{
  \item{data}{ Sample to be tested.}
  \item{mod0}{ The maximum number of modes in the null hypothesis. Default \code{mod0=1} (unimodality vs. multimodality test).}
  \item{method}{ The method employed for testing the number of modes. Available methods are: SI (Silverman, 1981), HY (Hall and York, 2001), FM (Fisher and Marron, 2001), HH (Hartigan and Hartigan, 1985), CH (Cheng and Hall,1998), NP (Ameijeiras--Alonso \emph{et al.}, 2016). Default \code{method="NP"}.}
  \item{B}{ Number of replicates used in the test. Default \code{B=500}.}
  \item{full.result}{ Logical. If \code{TRUE}, returns both the test statistic and the p-value. If \code{FALSE}, just the p-value is returned. Default is \code{FALSE}.
}
  \item{lowsup}{ Lower limit for the random variable support in the computation of the  critical bandwidth. Default is \code{-Inf}.}
  \item{uppsup}{ Upper limit for the random variable support in the computation of the critical bandwidth. Default is \code{Inf}.}
  \item{n}{ The number of equally spaced points at which the density is to be estimated. When n > 512, it is rounded up to a power of 2 as for \code{density} function. Default is \code{n=2^10} when \code{method} is SI, HY or FM and \code{n=2^15} when the \code{method} is CH or NP. Argument not used for other methods.}
  \item{tol}{ Accuracy for computing the critical bandwidth. Default \code{tol=10^(-5)}.}
  \item{tol2}{ Accuracy for integration of the calibration function in the \code{method} NP when the support is known. Default \code{tol2=10^(-5)}.}
  \item{methodhy}{ Different approaches when using Hall and York (2001) method. Default \code{methodhy=1}.}
  \item{alpha}{ Significance level employed for testing unimodality when method 1 of Hall and York (2001) is used. Default \code{alpha=0.05}.}
  \item{nMC}{ Number of Monte Carlo replicates used to approximate the p-value in the method 2 of Hall and York (2001). Default \code{nMC=100}.}
  \item{BMC}{ Number of bootstrap replicas used for computing the p--value in each Monte Carlo replicate of the Hall and York (2001) method 2. Default \code{BMC=100}.}
  }
\details{
The number of modes for the underlying density of a sample given by \code{data} can be tested with \code{modetest}. The null hypothesis states that the sample has \code{mod0} modes, and the alternative hypothesis is if it has more modes. The test used for calculating the p-value is specified in \code{method}. All the available proposals require bootstrap or Monte Carlo resamples (number specified in \code{B}).

Except when the support is employed, the typical usages are

\preformatted{
modetest(data,mod0=1,...)
modetest(data,mod0=1,method="NP",B=500,full.result=FALSE,n=NULL,tol=NULL)
}

Since a dichotomy algorithm is employed for computing the critical bandwidth (methods SI, HY, FM, NP), the parameter \code{tol} is used to determine a stopping time in such a way that the error committed in the computation of the critical bandwidth is less than \code{tol}.

The sample data can be perturbed in the methods using the excess of mass or the dip statistic (HH, CH and NP) in order to avoid important effects on the computation of the test statistic. See \code{excessmass}.

If a compact support containing a mode is known, it can be used to compute the Hall and York (2001) critical bandwidth. Note that in the case of the test proposed by Hall and York (2001), this support must be known. For their proposal, two methods are implemented. Method 1 is an asymptotic correction of Silverman (1981) test based on the limiting distribution of the test statistic. When \code{methodhy=1}, the significance level must be previously determined with \code{alpha}. Method 2 is based on Monte Carlo techniques. For this reason, when \code{methodhy=2}, the number of replicates \code{nMC} and the number of bootstrap replicas used for computing the p--value in each Monte Carlo replicate \code{BMC} are needed.

Typical usages are

\preformatted{
modetest(data,method="HY",lowsup=-1.5,uppsup=1.5,methodhy=1,alpha=NULL,...)
modetest(data,method="HY",B=500,full.result=FALSE,lowsup=-1.5,
uppsup=1.5,n=NULL,tol=NULL,methodhy=1,alpha=NULL)

modetest(data,method="HY",lowsup=-1.5,uppsup=1.5,methodhy=2,nMC=NULL,BMC=NULL,...)
modetest(data,method="HY",B=500,full.result=FALSE,lowsup=-1.5,
uppsup=1.5,n=NULL,tol=NULL,methodhy=2,nMC=NULL,BMC=NULL)
}

A modification of the proposal of Ameijeiras--Alonso \emph{et al.} (2016) can be also applied, by setting \code{method=NP} and including a known compact support for detecting the modes. The parameter \code{tol2} is the accuracy required in the integration of the calibration function. For more information, see Ameijeiras--Alonso \emph{et al.} (2016), the default approach, when the support is unknown, is given in Section 2.3 and, when it is provided, the approach shown in the Appendix B is employed.

Typical usages are

\preformatted{
modetest(data,mod0=1,lowsup=-1.5,uppsup=1.5,...)
modetest(data,mod0=1,B=500,full.result=FALSE,lowsup=-1.5,uppsup=1.5,
n=NULL,tol=NULL,tol2=NULL)
}


The NAs will be automatically removed.
}


\value{
  \item{p.value}{ P-value obtained after applying the test.}
  \item{statistic}{ Value of the test statistic. Critical bandwidth if the \code{method} is SI or HY; Cramer-von Mises statistic if the \code{method} is FM; the dip statistic if the \code{method} is HH, and the excess of mass when the \code{method} is CH or NP.}
}
\references{

Jose Ameijeiras--Alonso, Rosa M. Crujeiras, Alberto Rodríguez--Casal (2016). Mode testing, critical bandwidth and excess mass, \emph{arXiv preprint: 1609.05188}.

Cheng, M. Y. and Hall, P. (1998). Calibrating the excess mass and dip tests of modality, \emph{Journal of the Royal Statistical Society. Series B}, \bold{60}, 579--589.

Fisher, N.I. and Marron, J. S. (2001). Mode testing via the excess mass estimate,
 \emph{Biometrika}, \bold{88}, 419--517.

Hall, P. and York, M. (2001).  On the calibration of Silverman's test for multimodality, \emph{Statistica Sinica}, \bold{11}, 515--536.

Hartigan, J. A. and Hartigan, P. M. (1985). The Dip Test of Unimodality,
 \emph{Journal of the American Statistical Association}, \bold{86}, 738--746.

Silverman, B. W. (1981). Using kernel density estimates to investigate multimodality, \emph{Journal of the Royal Statistical Society. Series B}, \bold{43}, 97--99.

}

\author{ Jose Ameijeiras--Alonso, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

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



