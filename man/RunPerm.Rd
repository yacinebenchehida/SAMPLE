% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunPerm.R
\name{RunPerm}
\alias{RunPerm}
\title{Run permutations}
\usage{
RunPerm(input, replicates = 50)
}
\arguments{
\item{input}{Input dataframe (a dataframe object).}

\item{replicates}{Number of permutation replicates to perform (an integer; default replicates=50).}
}
\value{
A dataframe with the output values of the permutation analysis.
}
\description{
This function defines the number of permutations (default = 50) and subsequently runs the permutations on an input dataframe. A minimum number of 10 permutations is recommended in order to avoid having very large or very small stability values due to the stochastic nature of the process. On a large dataset, increasing the number of permutations can considerably slow down the analysis. Please refer to the Simulations section to see how changing permutation values can impact stability values.
}
\examples{
data("coral_symbionts")
perm <- RunPerm(input = coral_symbionts,replicates = 50)
perm

}
