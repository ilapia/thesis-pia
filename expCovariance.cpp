#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix expCovariance(NumericMatrix X1, NumericMatrix X2, double l, double sigma2, int trainingC) {
  // C = expCovariance(x1,x2,l,sigma2,training)
  // 
  // A function that takes matrices x1 and x2 of input locations, lengthscale
  // parameter l and magnitude parameter sigma2 and constructs a covariance
  // matrix between f(x1) and f(x2) corresponding to a Gaussian covariance
  // function. 
  //
  // Parameters are
  //  X1 : matrix of coordinates, n1 x d
  //  X2 : matrix of coordinates, n2 x d
  //  l : lengthscale
  //  sigma2 : variance parameter
  //  trainingC : an indicator telling whether we calculate the 
  //              training covariance matrix (1), that is X1 == X2
  //              or a generic matrix of covariances (0), that is X1 != X2
  
  // test that 
  if (X1.ncol() != X2.ncol()){
    throw std::range_error("number of columns in X1 and X2 do not match");
  }
  
  NumericMatrix C = NumericMatrix(X1.nrow(),X2.nrow());
  double d;
  if (trainingC == 1){
    // symmetric covariance matrix (training covariance)
    for(int ii = 0; ii < C.nrow(); ii++){ //C.nrow()-1 (upper triangular matrix without diagonal)
      for(int jj = ii+1; jj < C.ncol(); jj++){
        d = sqrt(sum(pow(X1.row(ii)-X2.row(jj),2))/pow(l,2));
        C(ii,jj) = sigma2*exp(-d);
        C(jj,ii) = sigma2*exp(-d);
      }
      C(ii,ii) = sigma2;
    }
  }
  else{
    // nonsymmetric matrix of covariances
    for(int ii = 0; ii < C.nrow(); ii++){
      for(int jj = 0; jj < C.ncol(); jj++){
        d = sqrt(sum(pow(X1.row(ii)-X2.row(jj),2))/pow(l,2));
        C(ii,jj) = sigma2*exp(-d);
      }
    }
  }
  return C;
}

