\encoding{utf8}
\name{multimode-package}
\alias{multimode-package}
\alias{multimode}
\docType{package}
\title{
Mode testing and exploring.
}

\description{
Different examples and methods for testing (including different proposals described in Ameijeiras-Alonso et al., 2016 <arXiv:1609.05188>) and exploring (including the mode tree, mode forest and SiZer map) the number of modes using nonparametric techniques.
}
\details{
\tabular{ll}{
Package: \tab multimode\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2018-06-08\cr
License: \tab GPL-3\cr
NeedsCompilation: \tab yes\cr
LazyData: \tab yes\cr
}
This package incorporates the function \code{modetest} which tests if the number of modes of a sample is equal to a given number (against if it is greater).
Functions \code{bw.crit} and \code{excessmass} provide the critical bandwidth and the excess mass statistic, respectively.
Function \code{nmodes} computes the number of modes for a given bandwidth. Given a certain number of modes, function \code{locmodes} provides the estimation of the locations of modes and antimodes and their density value.
Functions \code{modetree} and \code{modeforest} provide the mode tree and forest, respectively; they represent the estimated mode locations for different bandwidths. Function \code{sizer} can be used for determining where the smoothed curve is significantly increasing or decreasing. Registries with missing data are removed.
For a complete list of functions, use \code{library(help="multimode")}.
}

\references{
Ameijeiras--Alonso, J., Crujeiras, R.M. and Rodríguez--Casal, A. (2016). Mode testing, critical bandwidth and excess mass, \emph{arXiv preprint: 1609.05188}.

Ameijeiras--Alonso, J., Crujeiras, R.M. and Rodríguez--Casal, A. (2018). multimode: An R Package for Mode Assessment, \emph{arXiv preprint: 1803.00472}.

Chaudhuri, P. and Marron, J. S. (1999). SiZer for exploration of structures in curves, \emph{ Journal of the American Statistical Association}, \bold{94}, 807--823.

Minnotte, M. C., Marchette, D. J. and Wegman, E. J. (1998). The bumpy road to the mode forest, \emph{Journal of Computational and Graphical Statistics}, \bold{7}, 239--251.
}


%\keyword{ mode testing}
%\keyword{ kernel density estimation}
%\keyword{ mode estimation}

\section{Acknowledgements}{ This work has been supported by Projects MTM2016--76969--P (Spanish State Research Agency, AEI) and MTM2013--41383--P (Spanish Ministry of Economy, Industry and Competitiveness), both co--funded by the European Regional Development Fund (ERDF), IAP network (Developing crucial Statistical methods for Understanding major complex Dynamic Systems in natural, biomedical and social sciences, StUDyS) from Belgian Science Policy. Work of J. Ameijeiras-Alonso has been supported by the PhD grant BES-2014-071006 from the Spanish Ministry of Economy and Competitiveness.}




