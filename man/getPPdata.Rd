% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getPPdata.R
\name{getPPdata}
\alias{getPPdata}
\title{Data preprocessing}
\usage{
getPPdata(x, nfeatures = 2000)
}
\arguments{
\item{x}{an object of class \code{matrix}. Each column corresponds to a sample and each row to a variable.}

\item{nfeatures}{a numeric value represents the number of the filtered features.}
}
\value{
an object of class \code{matrix}.
}
\description{
This function is used to preprocess input data.
}
