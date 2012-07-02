\name{kaps}
\alias{kaps}
\title{K-adaptive partitioing for survival data}
\description{
    Conduct K-adaptive partitioning algorithm for censored survival data, especially with an application to SEER (Surveillance Epidemiology and End Results).
}
\usage{
	kaps(formula, data, K = 2:5, V = 5, mindat, ...)
}
\arguments{
  \item{formula}{an object with the type of Formula, with the response of the left hand side of the '~' operator, and the covariate terms on the right side. The response have be a survival object with survival time and censoring status by \code{Surv}. See \link[=Formula]{Formula}. }
  \item{data}{a dataset with variables used in \code{Formula}. It needs at least three variables including survival time, censoring status, and a covariate. The structure of this object is data.frame.}
  \item{K}{
    an object with the type of numeric to determine the number of subgroups.
  }
  \item{V}{
    an object with the type of numeric to determine the parameter V used in V-fold cross validation (CV). A default value is 5 (V=5). 
  }
  \item{mindat}{
    an object with the type of numeric in order to set the minimum number of observations at each subgroup.
  }
  \item{\ldots}{
	an object with the type of S4 class, "kapsOptions", that provides several minors paramters. See \link{kaps.control}.
  }
}
\details{
  This function provides routines for K-adaptive partitioning algorithm.
}
\author{
    Soo-Heang Eo \email{hanansh@korea.ac.kr} \cr
    HyungJun Cho \email{hj4cho@korea.ac.kr} \cr
}
\references{
Eo, SH, Hong, S-M and Cho, H (2012+). K-Adaptive Partitioning for Survival data with the application to SEER, submitted to Journal of Statitical Software.  
}
\seealso{
	 \code{\link{show}}, \code{\link{plot}}, \code{\link{predict}}, \code{\link{print}} and \code{\link{summary}} for the convenient use of kaps() \cr
	 \code{\link{kaps.control}} to control kaps() more detail \cr
	 \code{\link{lrtree}}, \code{\link{lrtree.control}} for counterpart of kaps()
}
\examples{
data(toy)
f <- Surv(time, status) ~ meta
# Fit maps algorithm without cross validation.
# It means the step to finding optimal K is not entered.
fit1 <- kaps(f, data = toy, K = 3)

# show the object of kaps (it contains apss S4 class)
fit1

# plot Kaplan-Meire survival curves
plot(fit1)

# To conduct cross-validation, use V argument without fold input.
fit2 <- kaps(f, data = toy, K= 2:4, V=5) # It conduct kaps with V-fold cross-validation when V=5.
fit2

plot(fit2)
print(fit2,K=2)
summary(fit2)
summary(fit2,K=2)

## compare with other algorithms
#fit3 <- lrtree(f, data = toy, minbucket = 10)
#fit3
#par(mfrow = c(1,1))
#plot(fit3)

# require(party)
# fit4 <- ctree(f, data = toy)
}
\keyword{kaps}
