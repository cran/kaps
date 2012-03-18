\name{lrtree}
\alias{lrtree}
\title{Recursive partitioning for censored survival data by log-rank test}
\description{
    Conduct recursive partitioning algorithm for censored survival data by log-rank related tests.
}
\usage{
	lrtree(Formula, data, subset = NULL, weights = NULL, ...) 
}

\arguments{
  \item{Formula}{an object with the type of Formula, with the response of the left hand side of the '~' operator, and the covariate terms on the right side. The response have be a survival object with survival time and censoring status by \code{Surv}. See \link[=Formula]{Formula}. }
  \item{data}{a dataset with variables used in \code{Formula}. It needs at least three variables including survival time, censoring status, and a covariate. The structure of this object is data.frame.}
  \item{subset}{
    an optional vector to be used as a subset of observations.
  }
  \item{weights}{
    an optional vector to be used as a weight parameter.
  }
  \item{\ldots}{
    a S4 class object specifying minior parameters for lrtree() using \link{lrtree.control} 
  }
}
\author{
    Soo-Heang Eo \email{hanansh@korea.ac.kr} \cr
}
\seealso{
	\code{\link{km.curve}}, \code{\link{labels}}, \code{\link{plot}}, \code{\link{show}}, \code{\link{summary}},  \code{\link{text}} for the convenient use of \code{lrtree} \cr	\code{\link{lrtree.control}} to control \code{lrtree} more detail \cr
	\code{\link{kaps}}, \code{\link{apss.control}} for counterpart of the lrtree
}
\examples{
## toy example
data(toy)
f <- Surv(time, status) ~ meta
fit <- lrtree(f, data = toy, B = 20)
### Print lrtree object
fit
### Print survival median time, 1 years, 3 years, and 5 years times. 
summary(fit)

### Plot tree diagram
plot(fit) # type = "naive"

### Plot Kaplan-Meire survival curves for terminal nodes.
km.curve(fit)

#f <- Surv(time, status) ~ ratio
#fit <- lrtree(f, data = toy, minbucket = 100, pre.pt = seq(0, 0.95, by = .05))
#plot(fit, type = "naive") # type = "naive"
### Print survival median time, i years, 3 years, and 5 years times by summary(). 
#summary(fit)

### Draw Kaplan-Meier survival curves
#km.curve(fit)
#legend(120,0.6, col=c(1:6)[fit.order],  lwd=2, 
#		   c("Node 8  (0.00 <= Ratio <= 0.10)",
#			 "Node 11 (0.45 <  Ratio <= 0.75)",
#			 "Node 3  (0.75 <  Ratio <= 1.00)"),
#		   bty="n")
}
\keyword{tree}
