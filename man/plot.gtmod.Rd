\name{plot.gtmod}
\alias{plot.gtmod}
\alias{print.gtmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot and print a \code{"gtmod"} object
}
\description{
Plot and print methods for \code{"gtmod"} objects, i.e., the results of the graphical tools for exploring the number of modes.
}
\usage{
\method{plot}{gtmod}(x, addplot = FALSE, xlab = NULL, ylab = NULL,
col.lines = "black", col.sizer = NULL, addlegend = TRUE,
poslegend = "topright", ...)

\method{print}{gtmod}(x,digits=getOption("digits"), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{ An object inheriting from class \code{"gtmod"}.}
    \item{addplot}{ Logical, if \code{TRUE} the \code{\link{modetree}} plot is added to the current one. Default is \code{FALSE}.}
    \item{xlab}{ A title for the x axis. See \code{title}.}
    \item{ylab}{ A title for the y axis. See \code{title}.}
    \item{col.lines}{ Colors employed in the \code{\link{modetree}}. If the length is two, the first color is employed for representing the location of the modes and the second one for the splitting of the modes. Default is \code{col.lines="black"}.}
    \item{col.sizer}{ Colors employed in the \code{\link{sizer}} for indicating the behaviour of the smoothed curve. The first color indicates where it is significantly increasing, the second where it is not significantly different from zero, the third where it is significantly negative and the forth where the data are too sparse for meaningful inference. Default \code{col.sizer=c("red","orchid","blue","grey")}.}
    \item{addlegend}{ Logical, if \code{TRUE} the legend is displayed in the \code{\link{sizer}}. Default \code{TRUE}.}
    \item{poslegend}{ Position where the legend in the \code{\link{sizer}} should be displayed. Default \code{posLegend="topright"}.}
    \item{digits}{ Number of significant digits to use, see \code{\link[base]{formatC}}.}
  \item{\dots}{ Arguments to be passed to subsequent methods, \code{\link[graphics]{plot.default}} (\code{\link{modetree}}) or \code{\link[graphics]{image.default}} (\code{\link{modeforest}} and \code{\link{sizer}}) for the \code{plot} method and \code{\link[base]{formatC}} for the \code{print} method.}
}

\author{ Jose Ameijeiras--Alonso, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\seealso{
\code{\link{modetree}}, \code{\link{modeforest}} and \code{\link{sizer}}, also for examples.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphical tools }% use one of  RShowDoc("KEYWORDS")
