\name{HyperParaControl-class}
\Rdversion{1.1}
\docType{class}
\alias{HyperParaControl-class}
\title{S4 Class \code{"HyperParaControl"}}
\description{
	S4 class of the object "HyperParaControl"
}
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("HyperParaControl", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{VarCtrl}:}{Object of class \code{"VarCtrl"} in which the variety of variable controls is set. }
    \item{\code{SplitCtrl}:}{Object of class \code{"SplitCtrl"} in which the variety of split controls is set. }
    \item{\code{TGCtrl}:}{Object of class \code{"TGCtrl"} in which the varieyt of tree growing controls is set. }

	}
}

\author{
Soo-Heang Eo <hanansh@korea.ac.kr>
}
\seealso{
\code{\link{lrtree.control}}
}
\examples{
	showClass("HyperParaControl")
}
\keyword{classes}