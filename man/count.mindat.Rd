\name{count.mindat}
\alias{count.mindat}
\title{Caculate the minimum sample size when the number of subgroups is given}
\description{
   This function caculates the minimum sample size of each partition when the number of subgroups is given. 
}
\usage{
	count.mindat(formula, data, part = 10)
}

\arguments{
  \item{formula}{a Formula object with a response on the left hand side of the '~' operator, and the covariate terms on the right side. The response has to be a survival object with survival time and censoring status in the \link[survival:Surv]{Surv} function. For more details, see \link[Formula:Formula]{Formula} page. }
  \item{data}{a data frame with variables used in formula. It needs at least three variables including survival time, censoring status, and a covariate. Multivariate covariates can be supported with "+" sign.}
  \item{part}{a numeric object to determine the number of subgroups we want to split.}
}
\seealso{
	\code{\link{kaps}}
}
\examples{
	data(toy)
	count.mindat(Surv(time,staus) ~ meta, data = toy, part =5)
	count.mindat(Surv(time,staus) ~ meta, data = toy, part =10)
	count.mindat(Surv(time,staus) ~ meta, data = toy)
}
\keyword{kaps}
