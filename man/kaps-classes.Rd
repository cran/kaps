\name{kaps-class}
\Rdversion{1.1}
\docType{class}
\alias{kaps-class}
\title{ Class \code{"kaps"}}
\description{
	A S4 class for \emph{K}-adaptive partitioning for survival data (kaps).
}
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("kaps")}. The most important slot is \code{groupID}, which is a vector consisting of the information about classified subgroups.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{evaluated function call}
    \item{\code{formula}:}{formula used in the model fitting}
    \item{\code{data}:}{data used in the model fitting}
    \item{\code{groupID}:}{information about the classified subgroup}
    \item{\code{index}:}{index for the optimal subgroup among the candidate K}
    \item{\code{X}:}{test statistic with the worst pair of subgroups for the split set}
    \item{\code{Z}:}{overall test statistic with K subgroups using the split set}
    \item{\code{pair}:}{selected pair of subgroups}
    \item{\code{split.var}:}{selected covariate in the model fitting}
    \item{\code{split.pt}:}{selected set of cut-off points}
    \item{\code{mindat}:}{minimum number of observations at a subgroup}
    \item{\code{test.stat}:}{Bonferroni corrected p-value matrix. The first row means overall p-values and the second one denotes p-values of the worst-pair against K. The column in the matrix describes the order of K.}
    \item{\code{over.stat.sample}:}{adjusted overall test statistic by Bootstrapping}
    \item{\code{pair.stat.sample}:}{adjusted worst-pair test statistic by Bootstrapping}
    \item{\code{groups}:}{candidate K used in the argument}
    \item{\code{results}:}{a list of results about each K}
    \item{\code{Options}:}{tuning parameters}
	}
}
\section{Methods}{
  \describe{
	\item{show}{\code{signature(object = "kaps")}: Same as the show method without the optional arguments }
	\item{print}{\code{signature(x = "kaps", K)}: Same as the print method with the specified number of subgroups K.}
	\item{plot}{\code{signature(x = "kaps", K)}: Plot an object}
	\item{predict}{\code{signature(object = "kaps")}: Predict an object by the estimated cut-off points}
	\item{summary}{\code{signature(object = "kaps")}: Summarize an object by survival times for each subgroup }
	}
}
\examples{
showClass("kaps")
}
\keyword{classes}
