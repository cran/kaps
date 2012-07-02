\name{lrtree}
\alias{lrtree}
\title{Recursive partitioning for censored survival data by log-rank test}
\description{
    Conduct recursive partitioning algorithm for censored survival data by log-rank related tests.
}
\usage{
	lrtree(Formula, data,  ...) 
}

\arguments{
  \item{Formula}{an object with the type of Formula, with the response of the left hand side of the '~' operator, and the covariate terms on the right side. The response have be a survival object with survival time and censoring status by \code{Surv}. See \link[=Formula]{Formula}. }
  \item{data}{a dataset with variables used in \code{Formula}. It needs at least three variables including survival time, censoring status, and a covariate. The structure of this object is data.frame.}
  \item{\ldots}{
    a S4 class object specifying minior parameters for lrtree() using \link{lrtree.control} 
  }
}
\author{
    Soo-Heang Eo \email{hanansh@korea.ac.kr} \cr
}
\seealso{
	\code{\link{km.curve}}, \code{\link{labels}}, \code{\link{plot}}, \code{\link{show}}, \code{\link{summary}},  \code{\link{text}} for the convenient use of \code{lrtree} \cr	\code{\link{lrtree.control}} to control \code{lrtree} more detail \cr
	\code{\link{kaps}}, \code{\link{kaps.control}} for counterpart of the lrtree
}
\examples{
## toy example
data(toy)
f <- Surv(time, status) ~ meta
fit <- lrtree(f, data = toy, eps = 1)
### Print lrtree object
fit
### Print survival median time, 1 years, 3 years, and 5 years times. 
summary(fit)

### Plot tree diagram
plot(fit) # type = "naive"
#plot(fit, type = "km") # survival curves for terminal nodes

### with multivariate covariates
fit1 <- lrtree(Surv(time, status) ~ meta + Race + Age + Grade, data = toy)
### determination of the tree size by cost-complexity pruning (0SE)
#fit2 <- prune(fit)
#fit2

### Plot Kaplan-Meire survival curves for terminal nodes.
#km.curve(fit1)
#legend(120,0.6, col=c(1:6),  lwd=2, 
#		   c("Node 2  (0 <= # of meta lymphs < 1)",
#			 "Node 3  (1 <= # of meta lymphs)"),
#		   bty="n")

# with Bootstrap stopping rule
#fit3 <- lrtree(f, data = toy, tree.size = "Boot")


}
\keyword{tree}
