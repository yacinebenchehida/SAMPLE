% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StabAssess.R
\name{StabAssess}
\alias{StabAssess}
\title{Assess the stability}
\usage{
StabAssess(data, stability_thres = 2, sucess_points = 50, diff = 1)
}
\arguments{
\item{data}{A dataframe generated by the function RunPerm()}

\item{stability_thres}{Threshold used to define stability (an integer). This value will be divided by the square root of the number of replicates (a float, default stability_thres = 2.0).}

\item{sucess_points}{Number of successive (mean) prevalence rates that are below a threshold (see parameter stability_thres) used to define stability (an integer; default stability_thres = 10).}

\item{diff}{Difference between the absolute minimum and maximum values among the all the means used to set the stability threshold (a float, default diff = 1.0).}
}
\value{
Two dataframes
}
\description{
Assess the stability
}
\examples{

data <- read.table("Example_file.txt")
perm = RunPerm(input = data,replicates = 50)
StabAssess(data = perm,stability_thres = 5 ,sucess_points = 5,diff = 2 )

}
