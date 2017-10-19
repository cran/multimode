\encoding{utf8}
\name{excessmass}
\alias{excessmass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Excess of mass }
\description{
   This function computes the excess of mass statistic.
}
\usage{
excessmass(data,mod0=1)
}
\arguments{
  \item{data}{ Sample for computing the excess of mass.}
  \item{mod0}{ Number of modes for which the excess of mass is calculated. Default \code{mod0=1}.}
  }
\details{
With \code{excessmass}, the excess of mass statistic, introduced by Müller and Sawitzki (1991), for the integer number of modes specified in \code{mod0} is computed.

If there are repeated data in the sample or the distance between different pairs of data points shows ties, a data perturbation is applied. This modification is made in order to avoid the discretization of the data which has important effects on the computation of the test statistic. The perturbed sample is obtained by adding a sample from the uniform distribution in minus/plus a half of the minimum of the positive distances between two sample points.

The NAs will be automatically removed.
}
\value{
  Excess of mass for \code{mod0} modes of the sample.
}
\references{

Müller, D. W. and Sawitzki, G. (1991). Excess mass estimates and tests for multimodality, \emph{The Annals of Statistics}, \bold{13}, 70--84.

}

\author{ Jose Ameijeiras--Alonso, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\examples{
# Excess of mass statistic for one mode
set.seed(2016)
data=rnorm(50)
excessmass(data)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mode testing }% __ONLY ONE__ keyword per line


