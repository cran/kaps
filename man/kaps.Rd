\name{kaps}
\alias{kaps}
\title{K-adaptive partitioing for survival data}
\description{
    Conduct \emph{K}-adaptive partitioning algorithm for survival data
}
\usage{
	kaps(formula, data, K = 2:4, mindat, type = c("perm", "NULL"), ...)
}
\arguments{
  \item{formula}{Formula object with a response on the left hand side of the '~' operator, and the covariate terms on the right side. The response has to be a survival object with survival time and censoring status in the \link[survival:Surv]{Surv} function. For more details, see \link[Formula:Formula]{Formula} page. }
  \item{data}{data frame with variables used in formula. It needs at least three variables including survival time, censoring status, and a covariate. Multivariate covariates can be supported with "+" sign.}
  \item{K}{number of subgroups used in the model fitting. The default value is 2:4 which means finding optimal subgroups ranging from 2 to 4. 
  }
  \item{type}{Select a type of algorithm in order to find optimal number of subgroups. Two options are provided: \code{perm} and \code{NULL}. The \code{perm} chooses subgroups using permutation procudures, while the \code{NULL} passes a optimal selection algorithm.
  }
  \item{mindat}{
    the minimum number of observations at each subgroup. The default value is 5\% of observations.
  }
  \item{\ldots}{
	a list of tuning parameters with the class, "kapsOptions". For more details, see \link{kaps.control}.
  }
}
\details{
  This function provides routines to conduct KAPS algorithm which is designed to classify cut-off values by the minimax-based rule. 
}
\value{
	The function returns an object with class "kaps" with the following slots.
    \item{\code{call}:}{evaluated function call}
    \item{\code{formula}:}{formula to be used in the model fitting}
    \item{\code{data}:}{data to be used in the model fitting}
    \item{\code{groupID}:}{information about the subgroup classification}
    \item{\code{index}:}{an index for the optimal subgroup among the candidate K}
    \item{\code{X}:}{test statistic with the worst pair of subgroups for the split set s}
    \item{\code{Z}:}{the overall test staitstic with K subgroups using the split set s}
    \item{\code{pair}:}{selected pair of subgroups}
    \item{\code{split.var}:}{selected covariate in the model fitting }
    \item{\code{split.pt}:}{selected set of cut-off points}
    \item{\code{mindat}:}{minimum number of observations at a subgroup}
    \item{\code{test.stat}:}{Bonferroni corrected p-value matrix. The first row means overall p-values and the second one denotes p-values of the worst-pair against K. The column in the matrix describes the order of K.}
    \item{\code{over.stat.sample}:}{adjusted overall test statistic by Bootstrapping}
    \item{\code{pair.stat.sample}:}{adjusted worst-pair test statistic by Bootstrapping}
    \item{\code{groups}:}{candidate K used in the argument}
    \item{\code{results}:}{list object about the results of each candidate K}
    \item{\code{Options}:}{tuning parameters}
}
\author{
    Soo-Heang Eo \email{eo.sooheang@gmail.com} \cr
    Seung-Mo Hong \email{smhong28@gmail.com} \cr
    HyungJun Cho \email{hj4cho@korea.ac.kr} \cr
}
\references{
S-H Eo, S-M Hong and H Cho (2014). K-adaptive partitioning for survival data, \emph{submitted}.
}
\seealso{
	 \code{\link{show}}, \code{\link{plot}}, \code{\link{predict}}, \code{\link{print}} and \code{\link{summary}} for the convenient use of kaps() \cr
	 \code{\link{kaps.control}} to control kaps() more detail \cr
   \code{\link{count.mindat}} to calculate minimum subgroup sample size 
}
\examples{
  \dontrun{
    data(toy)
    f <- Surv(time, status) ~ meta
    # Fit kaps algorithm without cross-validation.
    # It means the step to finding optimal K is not entered.
    fit1 <- kaps(f, data = toy, K = 3)

    # show the object of kaps (it contains apss S4 class)
    fit1

    # plot Kaplan-Meire estimates
    plot(fit1)

    # Fit kaps algorithm for selection optimal number of subgropus.
    fit2 <- kaps(f, data = toy, K= 2:4) 
    fit2

    # plot outputs with subgroup selection
    require(locfit) # for scatterplot smoothing
    plot(fit2)

    print(fit2,K=2)
    summary(fit2)
    summary(fit2,K=2)

    # require(party)
    # fit4 <- ctree(f, data = toy)
  }
}
\keyword{kaps}
