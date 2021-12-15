#' Generate simulation data (Complete data with group predictors)
#'
#' This function helps you quickly generate simulation data based on transformation model.
#' You just need to input the sample and dimension of the data
#' you want to generate and the covariance parameter pho.
#' This simulated example comes from Example 2 introduced by Li et al.(2012)
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
#' @importFrom stats rnorm
#' @export
#' @author Xuewei Cheng \email{xwcheng@csu.edu.cn}
#' @examples
#' n=100;
#' p=200;
#' pho=0.5;
#' data=GendataGP(n,p,pho)
#'
#' @references
#'
#' Li, R., W. Zhong, and L. Zhu (2012). Feature screening via distance correlation learning. Journal of the American Statistical Association 107(499), 1129–1139.
GendataGP <- function(n,p,rho)# n sample size; p dimension size.
{
  sig=matrix(0,p,p);
  sig=rho^abs(row(sig)-col(sig));
  diag(sig)<-rep(1,p);
  X=mvrnorm(n,rep(0,p),sig);
  q1=quantile(X[,3],0.25)
  q2=quantile(X[,3],0.5)
  q3=quantile(X[,3],0.75)
  Y=2*X[,1]+2*X[,2]+2*X[,4]+2*X[,5]
  +1*as.numeric(X[,3]<q1)+2*as.numeric(X[,3]>q1 & X[,3]<q2 )
  +3*as.numeric(X[,3]>q3)+rnorm(n)
  return(list(X=X,Y=Y));
}





