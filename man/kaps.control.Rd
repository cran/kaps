\name{kaps.control}
\alias{kaps.control}
\title{Control minor parameters for "apss" object}
\description{
    It provides the parameters for kaps(), which are less in use for practical situations.}
\usage{
	kaps.control(pre.pt = list(), scope = list(), rho = 0, ncl = 1, 
    fold = TRUE, lower.limit = 0, upper.limit = 12, 
	shortcut = TRUE, splits = c("logrank","maxstat"),
    p.adjust.methods = c("none","holm", "hochberg", "hommel", "bonferroni", 
        "BH", "BY", "fdr"))
}
\arguments{
  \item{pre.pt}{a parameter with the type of list in order to treat pre-specified split candidates. Use the option as list(var = split points), i.e., x = 1:100}
  \item{scope}{a parameter with the type of list in order to treat pre-determined split range. Use the option as list(var = ranges), i.e., x = c(1,100)}
  \item{rho}{a scalar parameter that controls the type of logrank test. See \link[=survival]{survdiff}.}
  \item{ncl}{a scalar parameter with the type of integer to determine the number of cores used in parallel computing}
  \item{fold}{a logical parameter to determine the use of V-fold Cross-Validation in the step of finding optimal K subgroups}
  \item{lower.limit}{a scalar parameter with the type of numeric to treat pre-determined overall lower bound. Default is 0.}
  \item{upper.limit}{a scalar parameter with the type of numeric to treat pre-determined overall upper bound. Default is 12.}
  \item{shortcut}{a logical value. If TRUE, we skip the off-diagonal matrix in pairwise-comparison to reduce computational cost. Default is TRUE.}
  \item{splits}{a character value to determine the kind of pairwise test. Default is logrank test.}
  \item{p.adjust.methods}{a character parameter to select the criteria for the multiple comparison.}
  }
\seealso{
	\code{\link{kaps}}
}
\keyword{kaps}
