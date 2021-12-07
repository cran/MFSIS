#' Generate simulation data (Discrete response data)
#'
#' This function helps you quickly generate simulation data based on linear model.
#' You just need to input the sample and dimension of the data
#' you want to generate and the covariance parameter pho.
#'
#' @param n Number of subjects in the dataset to be simulated. It will also equal to the
#' number of rows in the dataset to be simulated, because it is assumed that each
#' row represents a different independent and identically distributed subject.
#' @param p Number of predictor variables (covariates) in the simulated dataset.
#' These covariates will be the features screened by model-free procedures.
#' @param rho The correlation between adjacent covariates in the simulated matrix X.
#' The within-subject covariance matrix of X is assumed to has the same form as an
#' AR(1) auto-regressive covariance matrix, although this is not meant to imply
#' that the X covariates for each subject are in fact a time series. Instead, it is just
#' used as an example of a parsimonious but nontrivial covariance structure. If
#' rho is left at the default of zero, the X covariates will be independent and the
#' simulation will run faster.
#'
#' @return the list of your simulation data
#' @import MASS
#' @importFrom MASS mvrnorm
#' @importFrom stats runif
#' @export
#' @author Xuewei Cheng \email{xwcheng@csu.edu.cn}
#' @examples
#' n=100;
#' p=200;
#' pho=0.5;
#' data=gendata4(n,p,pho)
gendata4 <- function(n,p,rho)# n sample size; p dimension size.
{
  sig=matrix(0,p,p);
  sig=rho^abs(row(sig)-col(sig));
  diag(sig)<-rep(1,p);
  X=mvrnorm(n,rep(0,p),sig);
  b=c(0.6,0.6,0.6,0.6,-0.9*sqrt(2))
  myrates=abs(X[,1:5]%*%b+rnorm(n))
  Y=as.integer(runif(n,myrates,myrates+5))
  return(list(X=X,Y=Y));
}





