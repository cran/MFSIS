#' Model-Free Feature Screening Based on the Projection Correlation
#'
#' A model-free screening method is based on the projection correlation
#' which measures the dependence between two random vectors.
#' This projection correlation based method does not require specifying a
#' regression model, and applies to data in the presence of heavy tails
#' and multivariate responses. It enjoys both sure screening and
#' rank consistency properties under weak assumptions.
#'
#' @param X The design matrix of dimensions n * p. Each row is an observation vector.
#' @param Y The response vector of dimension n * 1.
#' @param nsis Number of predictors recruited by PCSIS. The default is n/log(n).
#'
#' @return the labels of first nsis largest active set of all predictors
#'
#' @import reticulate
#' @importFrom reticulate py_config
#' @importFrom reticulate source_python
#'
#' @export
#' @author Xuewei Cheng \email{xwcheng@csu.edu.cn}
#' @examples
#'
#'have_numpy=reticulate::py_module_available("numpy")
#'if (have_numpy){
#' req_py()
#' library(MFSIS)
#' n=100;
#' p=200;
#' pho=0.5;
#' data=gendata1(n,p,pho)
#' data=cbind(data[[1]],data[[2]])
#' colnames(data)[1:ncol(data)]=c(paste0("X",1:(ncol(data)-1)),"Y")
#' data=as.matrix(data)
#' X=data[,1:(ncol(data)-1)];
#' Y=data[,ncol(data)];
#' A=PCSIS(X,Y,n/log(n));A
#'}else{
#'    print('You should have the Python testing environment!')
#'}
#'
#' @references
#'
#' Zhu, L., K. Xu, R. Li, and W. Zhong (2017). Projection correlation between two random vectors. Biometrika 104(4), 829–843.
#'
#' Liu, W., Y. Ke, J. Liu, and R. Li (2020). Model-free feature screening and FDR control with knockoff features. Journal of the American Statistical Association, 1–16.
PCSIS=function(X,Y,nsis=(dim(X)[1])/log(dim(X)[1])){
  if (dim(X)[1]!=length(Y)) {
    stop("X and Y should have same number of rows!")
  }
  if (missing(X)|missing(Y)) {
    stop("The data is missing!")
  }
  if (TRUE%in%(is.na(X)|is.na(Y)|is.na(nsis))) {
    stop("The input vector or matrix cannot have NA!")
  }
  if (inherits(Y,"Surv")) {
    stop("PCSIS can not implemented with object  of Surv")
  }
  req_py()
  n=dim(X)[1];
  p=dim(X)[2];
  reticulate::py_config()
  py_path=system.file("python", "PCSIS.py", package="MFSIS")
  reticulate::source_python(py_path,envir=globalenv())
  corr=c();
  get_arccos_1d=NULL;
  projection_corr_1d=NULL;
  A_y=get_arccos_1d(Y)
  for (j in 1:p){
    A_x=get_arccos_1d(X[,j])
    corr[j]=projection_corr_1d(A_x,A_y,n);
  }
  A=order(corr,decreasing=TRUE)
  return (A[1:nsis])
}





