\name{plot.gtmod}
\alias{plot.gtmod}
\alias{print.gtmod}
\alias{summary.gtmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot, print and summarize a \code{"gtmod"} object
}
\description{
Plot, print and summarize methods for \code{"gtmod"} objects, i.e., the results of the graphical tools for exploring the number of modes.
}
\usage{
\method{plot}{gtmod}(x, addplot = FALSE, xlab = NULL, ylab = NULL,
col.lines = "black", col.sizer = NULL, addlegend = TRUE,
poslegend = "topright", \dots)

\method{print}{gtmod}(x,digits=getOption("digits"), \dots)

\method{summary}{gtmod}(object, bandwidths=TRUE, digits = getOption("digits"), width=1, levelmf=NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{ An object inheriting from class \code{"gtmod"}.}
    \item{object}{ An object of class \code{"gtmod"} for which a summary is desired.}
    \item{addplot}{ Logical, if \code{TRUE} the \code{\link{modetree}} plot is added to the current one. Default is \code{FALSE}.}
    \item{xlab}{ A title for the x axis. See \code{title}.}
    \item{ylab}{ A title for the y axis. See \code{title}.}
    \item{col.lines}{ Colors employed in the \code{\link{modetree}}. If the length is two, the first color is employed for representing the location of the modes and the second one for the splitting of the modes. Default is \code{col.lines="black"}.}
    \item{col.sizer}{ Colors employed in the \code{\link{sizer}} for indicating the behaviour of the smoothed curve. The first color indicates where it is significantly increasing, the second where it is not significantly different from zero, the third where it is significantly negative and the forth where the data are too sparse for meaningful inference. Default \code{col.sizer=c("red","orchid","blue","grey")}.}
    \item{addlegend}{ Logical, if \code{TRUE} the legend is displayed in the \code{\link{sizer}}. Default \code{TRUE}.}
    \item{poslegend}{ Position where the legend in the \code{\link{sizer}} should be displayed. Default \code{posLegend="topright"}.}
    \item{digits}{ Number of significant digits to use, see \code{\link[base]{formatC}}.}
    \item{bandwidths}{ Logical, if \code{TRUE} the bandwidths range for which the modes are estimated is showed. Default \code{TRUE}.}
    \item{width}{ The total field width in each showed number, see \code{\link[base]{formatC}}.}
    \item{levelmf}{ In the \code{\link{modeforest}}, just the pixels with a percentage (of times that an estimated mode falls on them) greater than this value will be considered. Default is \code{levelmf=0.5}.}
  \item{\dots}{ Arguments to be passed to subsequent methods, \code{\link[graphics]{plot.default}} (\code{\link{modetree}}) or \code{\link[graphics]{image.default}} (\code{\link{modeforest}} and \code{\link{sizer}}) for the \code{plot} method and \code{\link[base]{formatC}} for the \code{print} and the \code{summary} method.}
}

\author{ Jose Ameijeiras--Alonso, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\seealso{
\code{\link{modetree}}, \code{\link{modeforest}} and \code{\link{sizer}}, also for examples.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphical tools }% use one of  RShowDoc("KEYWORDS")
