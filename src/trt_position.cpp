#include "RcppArmadillo.h"
using namespace Rcpp;

//' Get the starting and ending position in a vector
//' 
//' @description
//' Function \code{trt_position()} gives the starting and ending index of 
//' covariate values in the \eqn{i_1^{th}} group of all the covariate values 
//' from all treatment groups that are together in a vector. E.g., covariate 
//' values in group 1 start from first value to the \eqn{n_1^{th}} value; those 
//' in group 2 start from \eqn{n_1 + 1} and end at the \eqn{(n_1 + n_2)^{th}} 
//' value. This function is for retrieving the position of an observation when 
//' the covariate values from all treatments are stored together in one vector.
//' 
//' This is a rewrite of the original \code{trt_position()} function from the
//' PSSIM library using Rcpp to decrease runtime. If you load both this library 
//' and PSSIM, load this one second so the faster functions will overwrite the 
//' slower ones.
//'
//' @param i1 An integer between 1 and \code{length(n)}.
//' @param n A vector of the sample sizes.
//' 
//' @return A 2-d vector, where the first gives the lower bound of treatment
//' \code{i1} and the second gives the upper bound.
//' 
//' @examples
//' i = 2; n = c(7, 8); trt_position(i, n)
//'  
//' @seealso \code{\link[=NPtest_indept]{NPtest_indept}} uses this function to 
//' run the independence test.
//' 
//' @references
//' Haiyan Wang, Siti Tolos, and Suojin Wang (2010). A Distribution Free
//'  Nonparametric Test to Detect Dependence Between a Response Variable and
//'  Covariate in Presence of Heteroscedastic Treatment Effects.
//'  The Canadian Journal of Statistics. 38(3), 408433. Doi:10.1002/cjs.10068
//' 
//' @useDynLib fastPSSIM
//' 
//' @export
// [[Rcpp::export]]
IntegerVector trt_position(int i1, IntegerVector n) {
  
  int lower;
  
  if (i1 == 1) {
    lower = 1;
  } else {
    lower = sum(n[Range(0, i1 - 2)]) + 1;
  } 
  int upper = sum(n[Range(0, i1 - 1)]);
  
  return IntegerVector::create(lower, upper);
} 