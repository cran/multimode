\encoding{utf8}
\name{excessmass}
\alias{excessmass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Excess mass }
\description{
   This function computes the excess mass statistic.
}
\usage{
excessmass(data,mod0=1,approximate=FALSE,gridsize=NULL,full.result=F)
}
\arguments{
  \item{data}{ Sample for computing the excess mass.}
  \item{mod0}{ Number of modes for which the excess mass is calculated. Default \code{mod0=1}.}
  \item{approximate}{ If this argument is TRUE then the excess mass value is approximated. Default \code{approximate=FALSE}.}
    \item{gridsize}{ When \code{approximate=TRUE}, number of endpoints at which the \eqn{C_m(\lambda)} sets are estimated (first element) and number of possible values of \eqn{\lambda} (second element). Default is \code{gridsize=c(20,20)}.}
  \item{full.result}{ If this argument is TRUE then it returns the full result list, see below. Default \code{full.result=FALSE}.}
}
\details{
With \code{excessmass}, the excess mass test statistic, introduced by Müller and Sawitzki (1991), for the integer number of modes specified in \code{mod0} is computed.

The excess mass test statistic for k modes is defined as \eqn{\max_{\lambda} \{D_{n,k+1}(\lambda)\}}, where \eqn{D_{n,k+1}(\lambda)=(E_{n,k+1}(P_n,\lambda)-E_{n,k}(P_n,\lambda))}. The empirical excess mass function for \eqn{k} modes is defined as \eqn{E_{n,k}(P_n,\lambda)=\sup_{C_1(\lambda),\ldots,C_k(\lambda)} \{\sum_{m=1}^k P_n (C_m(\lambda)) - \lambda ||C_m(\lambda)|| \}}, being the sets \eqn{C_m(\lambda)} closed intervals with endpoints the data points.

When \code{mod0>1} and the sample size is \emph{large}, a two--steps approximation (\code{approximate=TRUE}) can be performed in order to improve the computing time efficiency. First, since the possible \eqn{\lambda} candidates to maximize \eqn{D_{n,k+1}(\lambda)} can be directly obtained from the sets that maximize \eqn{E_{n,k+1}} and \eqn{E_{n,k}} (see Appendix E in Ameijeiras--Alonso \emph{et al.}, 2016), the possible values of \eqn{\lambda} are computed by looking to the empirical excess mass function in \code{gridsize[1]} endpoints candidates for \eqn{C_m(\lambda)} and also in the \eqn{\lambda} values associated to the empirical excess mass for one mode. Once a \eqn{\lambda} maximizing the approximated values of \eqn{D_{n,k+1}(\lambda)} is chosen, in order to obtain the approximation of the excess mass test statistic, in its neighborhood, a grid of possible values of \eqn{\lambda} is created, being its length equal to \code{gridsize[2]}, and the exact value of \eqn{D_{n,k+1}(\lambda)} is calculated for these values of \eqn{\lambda} (using the algorithm proposed by Müller and Sawitzki, 1991).


If there are repeated data in the sample or the distance between different pairs of data points shows ties, a data perturbation is applied. This modification is made in order to avoid the discretization of the data which has important effects on the computation of the test statistic. The perturbed sample is obtained by adding a sample from the uniform distribution in minus/plus a half of the minimum of the positive distances between two sample points.

The NAs will be automatically removed.
}
\value{
  Depending on \code{full.result} either a number, the excess mass statistic for \code{mod0} modes, or an object of class \code{"estmod"} which is a \code{list} containing the following components:
  \item{nmodes}{ The specified hypothesized value of the number of modes.}
  \item{sample.size}{ The number of non-missing observations in the sample used for computing the excess mass.}
  \item{excess.mass}{ Value of the excess mass test statistic.}
  \item{approximate}{ A logical value indicating if the excess mass was approximated.}
}
\references{

Jose Ameijeiras--Alonso, Rosa M. Crujeiras, Alberto Rodríguez--Casal (2016). Mode testing, critical bandwidth and excess mass, \emph{arXiv preprint: 1609.05188}.

Müller, D. W. and Sawitzki, G. (1991). Excess mass estimates and tests for multimodality, \emph{The Annals of Statistics}, \bold{13}, 70--84.

}

\author{ Jose Ameijeiras--Alonso, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\examples{
# Excess mass statistic for one mode
set.seed(2016)
data=rnorm(50)
excessmass(data)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mode testing }% __ONLY ONE__ keyword per line


