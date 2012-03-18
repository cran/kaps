\name{apss-class}
\Rdversion{1.1}
\docType{class}
\alias{apss-class}
\title{ S4 Class \code{"apss"}}
\description{
	S4 class of the object "apss"
}
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("apss", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{Object of class \code{"language"}  }
    \item{\code{formula}:}{Object of class \code{"Formula"}  }
    \item{\code{data}:}{Object of class \code{"data.frame"}  }
    \item{\code{where}:}{Object of class \code{"vector"}  }
    \item{\code{index}:}{Object of class \code{"integer"}  }
    \item{\code{Chisq}:}{Object of class \code{"numeric"}  }
    \item{\code{pvalue}:}{Object of class \code{"numeric"} }
    \item{\code{WH}:}{Object of class \code{"numeric"} }
    \item{\code{t}:}{Object of class \code{"numeric"} }
    \item{\code{Z}:}{Object of class \code{"numeric"} }
    \item{\code{pair}:}{Object of class \code{"numeric"} }
    \item{\code{split.var}:}{Object of class \code{"character"} }
    \item{\code{split.pt}:}{Object of class \code{"numeric"} }
    \item{\code{mindat}:}{Object of class \code{"numeric"} }
    \item{\code{elbow}:}{Object of class \code{"matrix"} }
    \item{\code{groups}:}{Object of class \code{"vector"} }
    \item{\code{candid}:}{Object of class \code{"list"} }
    \item{\code{Options}:}{Object of class \code{"apssOptions"} }
	}
}
\section{Methods}{
  \describe{
	\item{show}{\code{signature(object = "apss")}: Same as the show method without the optional arguments }
	\item{print}{\code{signature(object = "apss", K)}: Same as the print method with the specified number of subgroups K.}
	\item{plot}{\code{signature(object = "apss")}: Plots an apss object on the current graphics device}
	\item{predict}{\code{signature(object = "apss")}: Predict an apss object by the V-fold cross-validation}
	\item{summary}{\code{signature(object = "apss")}: Summarize an apss object with survival times for subgroup K }
	}
}
\author{
Soo-Heang Eo <hanansh@korea.ac.kr>
}
\examples{
showClass("apss")
}
\keyword{classes}
