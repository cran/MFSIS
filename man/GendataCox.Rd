% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GendataCox.R
\name{GendataCox}
\alias{GendataCox}
\title{Generate simulation data (Survival data based on the Cox model)}
\usage{
GendataCox(n, p, rho, beta = c(rep(1, 5), rep(0, p - 5)), lambda = 0.1)
}
\arguments{
\item{n}{Number of subjects in the dataset to be simulated. It will also equal to the
number of rows in the dataset to be simulated, because it is assumed that each
row represents a different independent and identically distributed subject.}

\item{p}{Number of predictor variables (covariates) in the simulated dataset.
These covariates will be the features screened by model-free procedures.}

\item{rho}{The correlation between adjacent covariates in the simulated matrix X.
The within-subject covariance matrix of X is assumed to has the same form as an
AR(1) auto-regressive covariance matrix, although this is not meant to imply
that the X covariates for each subject are in fact a time series. Instead, it is just
used as an example of a parsimonious but nontrivial covariance structure. If
rho is left at the default of zero, the X covariates will be independent and the
simulation will run faster.}

\item{beta}{A vector with length of n, which are the coefficients that you want to generate
about Linear model. The default is beta=(1,1,1,1,1,0,...,0)^T;}

\item{lambda}{This parameter control the censoring rate in survival data.
The censored time is generated by exponential distribution with mean 1/lambda. The default
is lambda=0.1.}
}
\value{
the list of your simulation data
}
\description{
This function helps you quickly generate simulation data based on the Cox model.
You just need to input the sample and dimension of the data
you want to generate and the covariance parameter rho.
}
\examples{
n <- 100
p <- 200
rho <- 0.5
data <- GendataCox(n, p, rho)

}
\references{
Cox DR (1972). “Regression models and life-tables.” Journal of the Royal Statistical Society:Series B (Methodological), 34(2), 187–202.
}
\author{
Xuewei Cheng \email{xwcheng@hunnu.edu.cn}
}
