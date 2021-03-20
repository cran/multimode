\encoding{utf8}
\name{sizer}
\alias{sizer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ SIgnificant ZERo crossing }
\description{
   This function provides the SiZer map.
}
\usage{
sizer(data,method=2,bws=NULL,gridsize=NULL,alpha=0.05,B=NULL,n0=NULL,
cbw1=NULL,cbw2=NULL,display=TRUE,logbw=TRUE,from=NULL,to=NULL,logbw.regulargrid=NULL,\dots)
}
\arguments{
  \item{data}{ Sample in which the SiZer map is computed.}
  \item{method}{ The method employed for computing the SiZer map. Available methods are: 1 (q1, pointwise Gaussian quantiles), 2 (q2, approximate simultaneous over x Gaussian quantiles), 3 (q3, bootstrap quantile simultaneous over x) and 4 (q4, bootstrap quantile simultaneous over x and h). Default \code{method=2}.}
  \item{bws}{ Vector or range of bandwidths. If it is a vector of size two, then it is used a grid of bandwidths between the given values. Default lower bandwidth is twice the grid size used for estimating the density and upper bandwidth equal to the range of the data. Unless it is specified a vector of size greater than two, the number of bandwidths employed is equal to the second element of \code{gridsize}.}
  \item{gridsize}{ Number of grid points in the horizontal (values of the variable, first element) and vertical (bandwidths, second element) axis. Default is \code{c(512,151)}.}
    \item{alpha}{ Significance level employed for determining the significant features. Default \code{alpha=0.05}.}
  \item{B}{ Number of replicates used for generating the SiZer map when method q3 or q4 are used. Default \code{B=100}.}
  \item{n0}{ When the effective sample size is below this quantity, the pixel in the SiZer map is shaded grey. Default \code{n0=5}.}
  \item{cbw1}{ Number of modes for which the first critical bandwidth is calculated. This is the first bandwidth used to compute the SiZer map when \code{bws} is not specified.}
  \item{cbw2}{ Number of modes for which the second critical bandwidth is calculated. This is the last bandwidth used to compute the SiZer map when \code{bws} is not specified.}
  \item{display}{ Logical, if \code{TRUE} the SiZer map is plotted. Default \code{TRUE}.}
  \item{logbw}{ Logical, if \code{TRUE} the plot displays and returns the log10 bandwidths. Default \code{logbw=TRUE}.}
    \item{from}{ First plotted value of the variable. Default is below the data minimum.}
    \item{to}{ Last plotted value of the variable. Default is above the data maximum.}
    \item{logbw.regulargrid}{ Logical, if \code{TRUE} a regular grid of bandwidths is created over the log10 scale. Default \code{logbw.regulargrid=FALSE}.}
    \item{\dots}{ Arguments to be passed to subsequent methods, \code{\link{plot.gtmod}} for the \code{plot}, \code{print} and \code{summary} methods.}
}

\details{

With this function the assessment of SIgnificant ZERo crossing of the derivative of the smoothed curve are computed for the sample given in \code{data}. For this calculation, a kernel density estimator with Gaussian kernel and bandwidths \code{bws} is used. When it is plotted, at a given location (horizontal axis) and using a specified bandwidth parameter (vertical axis), the SiZer map represents where the curve is significantly increasing (blue color by default), decreasing (red by default) or not significantly different from zero (orchid, a light tone of purple, by default). Thus, for a given bandwidth, a region significantly increasing followed by a region significantly decreasing (blue-red pattern by default) indicates where there is a significant peak.

For methods q2, q3 and q4, it is calculated where the data are too sparse for meaningful inference (grey color by default). A location-bandwidth pixel is classified in this last category when the estimated Effective Sample Size is less than \code{n0}. For more information, see Chaudhuri and Marron (1999).

For methods q3 and q4, the bootstrap quantiles are computed generating \code{B} random samples drawn with replacement from \code{data}.

The NAs will be automatically removed.
}
\value{
  A list with class \code{"gtmod"} containing the following components:
  \item{sizer}{ Matrix indicating the significant behavior of the smoothed curve in each location-bandwidth pixel. One indicates significantly decreasing; two, not significantly different from zero; three, significantly increasing and four where the data are too sparse for meaningful inference.}
  \item{method}{ A number indicating what type of quantile was performed.}
  \item{lower.CI}{ Matrix containig the lower limit of the confidence interval in each location-bandwidth pixel.}
    \item{estimate}{ Matrix containing the derivative values of the kernel density estimation in each location-bandwidth pixel.}
  \item{upper.CI}{ Matrix containig the upper limit of the confidence interval in each location-bandwidth pixel.}
  \item{ESS}{ Matrix containing the Effective Sample Size in each location-bandwidth pixel.}
  \item{range.x}{ Employed location values to represent the SiZer map.}
  \item{range.bws}{ Employed bandwidths to compute the SiZer map.}
    \item{logbw}{ Logical value indicating if the bandwidths are given in the log10 scale.}
  \item{sample.size}{ The number of non-missing observations in the sample used for obtaining the SiZer map.}
  \item{call}{ The unevaluated expression, which consists of the named function applied to the given arguments.}
}
\references{

Ameijeiras-Alonso, J., Crujeiras, R.M. and Rodríguez-Casal, A. (2021). multimode: An R Package for Mode Assessment, \emph{Journal of Statistical Software}, \bold{97}, 1--32.

Chaudhuri, P. and Marron, J. S. (1999). SiZer for exploration of structures in curves, \emph{ Journal of the American Statistical Association}, \bold{94}, 807--823.

}


\author{ Jose Ameijeiras-Alonso, Rosa M. Crujeiras and Alberto Rodríguez-Casal}

\examples{
#SiZer map using a grid of bandwidths between 1 and 10
data(geyser)
data=geyser
sizer(data,bws=c(1,10))

\dontrun{
#Different methods for calculating the confidence limits

#Pointwise Gaussian quantiles
sizer(data,method=1,bws=c(1,10))

#Approximate simultaneous over x Gaussian quantiles
sizer(data,method=2,bws=c(1,10))

#Bootstrap quantile simultaneous over x
sizer(data,method=3,bws=c(1,10))

#Bootstrap quantile simultaneous over x and h
sizer(data,method=4,bws=c(1,10))
}

#Adding the original mode tree for this sample
modetree(data,bws=c(0.8,10),logbw=TRUE,addplot=TRUE,col.lines="white")

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphical tools }% __ONLY ONE__ keyword per line

