% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability.R
\name{stability}
\alias{stability}
\title{Assess the stability}
\usage{
stability(data, stability_thresh = 2, success_points = 50, diff = 1)
}
\arguments{
\item{data}{A dataframe generated by the function RunPerm().}

\item{stability_thresh}{Threshold used to define stability (an integer). This value will be divided by the square root of the number of replicates (a float; default stability_thresh = 2.0).}

\item{success_points}{Number of successive (mean) prevalence rates that are below a threshold (see argument stability_thresh) used to define stability (an integer; default stability_thresh = 10).}

\item{diff}{Difference between the absolute minimum and maximum values among the all means used to set the stability threshold (a float; default diff = 1.0).}
}
\value{
A dataframe with the output values of the stability assessment analysis.
}
\description{
This function determines the stability values for all of the input variables in the dataset (e.g. species, habitats). Stability here is defined as the minimum occurrence value at which the variable in question stayed stable. See the arguments below to understand how they can affect the stability values.
}
\examples{

data("coral_symbionts")
set.seed(812)
perm <-  RunPerm(input = coral_symbionts,replicates = 50)
stable <- stability(data = perm,stability_thresh = 5 ,success_points = 5,diff = 2 )
stable

}
