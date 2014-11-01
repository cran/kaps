\name{km.curve}
\alias{km.curve}
\title{Plot Kaplan-Meire survival curves}
\description{
    Plot a Kaplan-Meire survival curve for terminal nodes or selected subgroups.
}
\usage{
	km.curve(object, 
	x.lab = c(0,24,48,72,96,120, 144, 168, 192, 216, 240), lwd = 1.5, ...)
}

\arguments{
  \item{object}{an object from \code{kaps}}
  \item{x.lab}{X labels specified as arguments}
  \item{lwd}{line width}
  \item{\dots}{other arguments for plot object. See \link[graphics:plot]{plot} for details.}
}

\details{
  This function provides Kaplan-Meire survival curves with the estimated subgroups by \code{\link{kaps}}. 
}
\seealso{
	\code{\link{kaps}}
}
\keyword{methods}
