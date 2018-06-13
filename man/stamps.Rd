\encoding{utf8}
\name{stamps}
\alias{stamps}
\alias{stamps1}
\alias{stamps2}
\alias{stampstable}
\docType{data}
\title{ Stamps thickness
}
\description{
This dataset, analysed in Izenman and Sommer (1988) and Ameijeiras--Alonso \emph{et al.} (2017), consists of thickness measurements (in millimeters) of 485 unwatermarked used white wove stamps of the 1872 Hidalgo stamp issue of Mexico. All of them had an overprint with the year (1872 or either an 1873 or 1874) and some of them were watermarked (Papel Sellado or LA+-F).
}
\usage{
data(stamps)
data(stamps1)
data(stamps2)
data(stampstable)
}
\format{
  \code{stamps} includes the thickness (in millimeters) of the different stamps in the 1872 Hidalgo stamp issue. \code{stampstable} reproduces the original frequency table in Izenman and Sommer (1988), including the \code{thickness}, the overprinted years (\code{1872} or \code{1873-1874}) and the watermarks (\code{Papel Sellado} or \code{LA+-F}). \code{stamps1} includes both the \code{thickness} and the overprinted \code{year}. \code{stamps2} includes both the \code{thickness} and the \code{watermark}. Note that the stamp with a thickness equal to 0.118 and the watermark \code{"Papel Sellado"} was removed in \code{stamps2} since it does not have a corresponding year, probably this watermark is associated with the observation of \code{thickness} 0.117 or 0.119.
}
\details{
This is a classic example for determining the number of modes.
}
\source{
Obtained from the original table in Izenman and Sommer (1988).
}
\references{

Ameijeiras--Alonso, J., Crujeiras, R.M. and Rodríguez--Casal, A. (2016). Mode testing, critical bandwidth and excess mass, \emph{arXiv preprint: 1609.05188}.

Izenman, A. J., and Sommer, C. J. (1988). Philatelic mixtures and multimodal densities. \emph{Journal of the American Statistical association}, \bold{83}, 941--953.

}
\examples{
data("stamps")
# Histogram of Wilson (bin width 0.008)
seqx=seq(0.0585,0.1385,by=0.008)
hist(stamps,breaks=seqx)
# Histogram of Figure 1, Izenman and Sommer (bin width 0.002)
seqx=seq(0.0585,0.1385,by=0.002)
hist(stamps,breaks=seqx)

# Kernel density estimation of Izenman and Sommer (number of modes: 7)
locmodes(stamps,mod0=7,lowsup=0.04,uppsup=0.15,display=TRUE)
# Kernel density estimation of Ameijeiras-Alonso et al. (number of modes: 4)
locmodes(stamps,mod0=4,lowsup=0.04,uppsup=0.15,display=TRUE)

}
\keyword{datasets}

