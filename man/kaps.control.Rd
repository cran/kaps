\name{kaps.control}
\alias{kaps.control}
\title{Control tuning parameters for "kaps" object}
\description{
    Tuning parameters for an object "kaps"}
\usage{
	kaps.control(pre.pt = list(), scope = list(), 
  rho = 0, V = 5, ncl = 1, lower.limit = 0, upper.limit = 100,
  shortcut = TRUE, N.perm = 9999, N.boot = 200, alpha = 0.05,
  splits = c("logrank", "exact"),
  correct = c("Adj.Bf", "Bf", "None"),
  p.adjust.methods = c("none", "holm", "hochberg", 
    "hommel", "bonferroni",  "BH", "BY", "fdr"))
}
\arguments{
  \item{pre.pt}{list parameter that treats pre-specified split candidates. Use the option as list(var = split points), i.e., x = 1:100}
  \item{scope}{list parameter that treats pre-determined split range. Use the option as list(var = ranges), i.e., x = c(1,100)}
  \item{rho}{scalar parameter that controls the type of logrank test. See \link[=survdiff]{survdiff}.}
  \item{V}{numeric parameter that determines the number of folds in the cross-validation subgroup selection.}
  \item{ncl}{integer parameter that determines the number of cores to improve computing power}
  \item{lower.limit}{numeric parameter that treats pre-determined overall lower bound. Default is 0.}
  \item{upper.limit}{numeric parameter that treats pre-determined overall upper bound. Default is 100.}
  \item{shortcut}{logical parameter. If shortcut = TRUE, we skip the off-diagonal matrix in pairwise-comparison to reduce computational cost. The default value is TRUE.}
  \item{N.perm}{numeric parameter that gives the number of permutation samples used in the kaps algorithm. The default value is 9999.}
  \item{N.boot}{numeric parameter that gives the number of Bootstrap samples used in the bootstrap and permuting kaps algorithm. The default value is 200.}
  \item{alpha}{numeric parameter that provides a significant level in the process of Bootstrap and permuting algorithm.}
  \item{splits}{character parameter that detemines the kind of pairwise test. Default is logrank test. At this stage, the option \code{exact} is not working.}
  %\item{sel}{a character value to determine the kind of selection methods for the pairwise test stestatic in Bootstrap algorithm. The \code{mean} and \code{median} denote the Bootstrap meand and median, \code{trim} describes censoring-related trimmed mean, and \code{test} is the method using proportional test for the Bootsrap samples. In \code{test}, the null hypothesis is p = 0.05.}
  \item{correct}{character parameter to select the criteria for the multiple comparison in the simple permuting kaps algorithm.}
  \item{p.adjust.methods}{character parameter to select the criteria for the multiple comparison.}
  }
\seealso{
	\code{\link{kaps}}
}
\keyword{kaps}
