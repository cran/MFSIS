#' The Kolmogorov filter for variable screening
#'
#' A new model-free screening method called the fused Kolmogorov filter is proposed for high-dimensional data analysis.
#' This new method is fully nonparametric and can work with many types of covariates and response variables, including continuous,
#' discrete and categorical variables.
#'
#' @param X The design matrix of dimensions n * p. Each row is an observation vector.
#' @param Y The response vector of dimension n * 1.
#' @param nsis Number of predictors recruited by SIS. The default is n/log(n).
#' @param response.type The type of response, you can choose "continuous", "discrete" or "categorical". The default is continuous.
#' @param method The method you want to employ in screening procedure, you can choose "single" or "fused".
#'
#' @return the labels of first nsis largest active set of all predictors
#' @export
#' @author Xuewei Cheng \email{xwcheng@csu.edu.cn}
#' @examples
#'
#'n=100;
#'p=200;
#'pho=0.5;
#'data=gendata1(n,p,pho)
#'data=cbind(data[[1]],data[[2]])
#'colnames(data)[1:ncol(data)]=c(paste0("X",1:(ncol(data)-1)),"Y")
#'data=as.matrix(data)
#'X=data[,1:(ncol(data)-1)];
#'Y=data[,ncol(data)];
#'A=Kfilter(X,Y,n/log(n));A
#'
#' @references
#'
#' Mai, Q., & Zou, H. (2013). The Kolmogorov filter for variable screening in high-dimensional binary classification. Biometrika, 100(1), 229-234.
#'
#' Mai, Q., & Zou, H. (2015). The fused Kolmogorov filter: A nonparametric model-free screening method. The Annals of Statistics, 43(4), 1471-1497.
Kfilter<-function(X,Y,nsis=(dim(X)[1])/log(dim(X)[1]),response.type="continuous",method="fused"){
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
    stop("Kfilter can not implemented with object  of Surv")
  }
  N=NULL;
  nslices=NULL;
  slicing.scheme=NULL;
  method=match.arg(arg=method,choices=c("fused","single"))
  response.type=match.arg(arg=response.type,choices=c("continuous","discrete","categorical"))

  if(!is.null(slicing.scheme)){
    if(method=="fused")warning("A slicing scheme is given. Using a single Kolmogorov filter.")
    obj<-Kfilter_single(X=X,Y=Y,nsis)
  }else{
    if(response.type=="categorical"){
      Y<-factor(Y)
      obj<-Kfilter_single(X=X,Y=Y,nsis)
    }
    if(response.type=="continuous"|response.type=="discrete"){
      n<-nrow(X)
      if(is.null(N)){N<-ceiling(log(n))-2}
      if(is.null(nslices)){nslices<-3:(N+2)}
      obj<-Kfilter_fused(X=X,Y=Y,nsis)
    }
  }
  return (obj)
}


#' The Kolmogorov filter for variable screening in high-dimensional binary classification
#'
#' @param X The design matrix of dimensions n * p. Each row is an observation vector.
#' @param Y The response vector of dimension n * 1.
#' @param nsis Number of predictors recruited by SIS. The default is n/log(n).
#'
#' @return the labels of first nsis largest active set of all predictors
#' @importFrom stats quantile
#' @importFrom stats ks.test
#' @export
#'
#' @examples
#' n=100;
#' p=200;
#' pho=0.5;
#' data=gendata3(n,p,pho)
#' data=cbind(data[[1]],data[[2]])
#' colnames(data)[1:ncol(data)]=c(paste0("X",1:(ncol(data)-1)),"Y")
#' data=as.matrix(data)
#' X=data[,1:(ncol(data)-1)];
#' Y=data[,ncol(data)];
#' A=Kfilter_single(X,Y,n/log(n));A
#'
#' @references
#' #' Mai, Q., & Zou, H. (2013). The Kolmogorov filter for variable screening in high-dimensional binary classification. Biometrika, 100(1), 229-234.
Kfilter_single<-function(X,Y,nsis=(dim(X)[1])/log(dim(X)[1])){
  Y=as.factor(Y)
  Y.dm<-Y
  Y.dm<-as.numeric(Y.dm)
  K<-length(unique(Y.dm))
  p<-ncol(X)
  ks.stat<-matrix(0,p,K*(K-1)/2)

  nclass<-0
  for(j in 1:(K-1)){
    for(l in (j+1):K){
      nclass<-nclass+1
      for(i in 1:p){
        ks.stat[i,nclass]<-ks.test(X[Y.dm==j,i],X[Y.dm==l,i])$statistic}
    }
  }
  ks.stat.max0<-apply(ks.stat,1,max)
  k.rank<-rank(-ks.stat.max0,ties.method="max")
  return (k.rank[1:nsis])
}


#' The fused kolmogorov filter: a nonparametric model-free screening method
#'
#' @param X The design matrix of dimensions n * p. Each row is an observation vector.
#' @param Y The response vector of dimension n * 1.
#' @param nsis Number of predictors recruited by SIS. The default is n/log(n).
#'
#' @return the labels of first nsis largest active set of all predictors
#' @importFrom stats quantile
#' @importFrom stats ks.test
#' @export
#'
#' @examples
#' ##Scenario 1  generate discrete response data
#' n=100;
#' p=200;
#' pho=0.5;
#' data=gendata4(n,p,pho)
#' data=cbind(data[[1]],data[[2]])
#' colnames(data)[1:ncol(data)]=c(paste0("X",1:(ncol(data)-1)),"Y")
#' data=as.matrix(data)
#' X=data[,1:(ncol(data)-1)];
#' Y=data[,ncol(data)];
#' A1=Kfilter_fused(X,Y,n/log(n));A1
#'
#' ##Scenario 2  generate continuous response data
#' n=100;
#' p=200;
#' pho=0.5;
#' data=gendata1(n,p,pho)
#' data=cbind(data[[1]],data[[2]])
#' colnames(data)[1:ncol(data)]=c(paste0("X",1:(ncol(data)-1)),"Y")
#' data=as.matrix(data)
#' X=data[,1:(ncol(data)-1)];
#' Y=data[,ncol(data)];
#' A2=Kfilter_fused(X,Y,n/log(n));A2
#'
#' @references
#'  Mai, Q., & Zou, H. (2015). The fused Kolmogorov filter: A nonparametric model-free screening method. The Annals of Statistics, 43(4), 1471-1497.
Kfilter_fused<-function(X,Y,nsis=(dim(X)[1])/log(dim(X)[1])){

  n<-nrow(X)
  p<-ncol(X)

  N<-ceiling(log(n))-2
  nslices<-3:(N+2)

  ks.stat.single<-matrix(0,N,p)
  ks.stat.max<-rep(0,p)

  for(K in nslices){
    slicing.scheme<-quantile(Y,seq(0,1,1/K))
    slicing.scheme[1]<-slicing.scheme[1]-1
    slicing.scheme[K+1]<-slicing.scheme[K+1]+1
    Y.dm<-cut(Y,slicing.scheme,labels=c(1:K),right=F)
    ks.stat<-matrix(0,p,K*(K-1)/2)

    nclass<-0
    for(j in 1:(K-1)){
      for(l in (j+1):K){
        nclass<-nclass+1
        for(i in 1:p){
          ks.stat[i,nclass]<-ks.test(X[Y.dm==j,i],X[Y.dm==l,i])$statistic}
      }
    }
    ks.stat.max0<-apply(ks.stat,1,max)
    ks.stat.single[K-2,]<-ks.stat.max0
    ks.stat.max<-ks.stat.max+ks.stat.max0
  }
  k.rank=rank(-ks.stat.max,ties.method="max")
  return (k.rank[1:nsis])
}






