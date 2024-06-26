% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WLS.R
\name{WLS}
\alias{WLS}
\title{A Model-free Variable Screening Method Based on Leverage Score}
\usage{
WLS(X, Y, nsis = (dim(X)[1])/log(dim(X)[1]))
}
\arguments{
\item{X}{The design matrix of dimensions n * p. Each row is an observation vector.}

\item{Y}{The response vector of dimension n * 1.}

\item{nsis}{Number of predictors recruited by WLS. The default is n/log(n).}
}
\value{
the labels of first nsis largest active set of all predictors.
}
\description{
An innovative and effective sampling scheme based on leverage scores via singular value decompositions
has been proposed to select rows of a design matrix as a surrogate of the full data in linear regression.
Analogously, variable screening can be viewed as selecting rows of the design matrix. However, effective
variable selection along this line of thinking remains elusive. This method propose a
weighted leverage variable screening method by using both the left and right singular vectors of the design matrix.
}
\examples{

n <- 100
p <- 200
rho <- 0.5
data <- GendataLM(n, p, rho, error = "gaussian")
data <- cbind(data[[1]], data[[2]])
colnames(data)[1:ncol(data)] <- c(paste0("X", 1:(ncol(data) - 1)), "Y")
data <- as.matrix(data)
X <- data[, 1:(ncol(data) - 1)]
Y <- data[, ncol(data)]
A <- WLS(X, Y, n / log(n))
A

}
\references{
Zhong, W., Liu, Y., & Zeng, P. (2021). A Model-free Variable Screening Method Based on Leverage Score. Journal of the American Statistical Association, (just-accepted), 1-36.
}
\author{
Xuewei Cheng \email{xwcheng@hunnu.edu.cn}
}
