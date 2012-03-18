\name{Tree-class}
\Rdversion{1.1}
\docType{class}
\alias{Tree-class}
\title{ S4 Class \code{"Tree"}}
\description{
	S4 class of the object "Tree"
}
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("Tree", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{Object of class \code{"language"}  }
    \item{\code{formula}:}{Object of class \code{"Formula"}  }
    \item{\code{data}:}{Object of class \code{"data.frame"}  }
    \item{\code{frame}:}{Object of class \code{"data.frame"}  }
    \item{\code{controls}:}{Object of class \code{"HyperParaControl"}  }
    \item{\code{where}:}{Object of class \code{"vector"}  }
    \item{\code{subset}:}{Object of class \code{"numeric"} }
    \item{\code{weights}:}{Object of class \code{"numeric"}  }
    \item{\code{alpha}:}{Object of class \code{"numeric"} }
    \item{\code{alpha.prime}:}{Object of class \code{"numeric"}  }
	}
}
\section{Methods}{
  \describe{
	\item{summary}{\code{signature(object = "Tree")}: summarize the "Tree" object with survival median time for each terminal node. }
	\item{labels}{\code{signature(object = "Tree")}: Create Split Labels for an "Tree". }
	\item{show}{\code{signature(object = "Tree")}: Same as the print method without the optional arguments. }
	\item{plot}{\code{signature(object = "Tree")}: Plots a Tree object on the current graphics device.}
	\item{text}{\code{signature(object = "Tree")}: Labels the current graphic device of the tree dendrogram with text. }
	}
}
\author{
Soo-Heang Eo <hanansh@korea.ac.kr>
}
\examples{
showClass("Tree")
}
\keyword{classes}
