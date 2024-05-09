#include "RcppArmadillo.h"
using namespace Rcpp;

//' Index in one vector mapped to treatment and observation index
//' 
//' @description
//' Function \code{mapindex()} maps the 1-d index \eqn{r = 1, \ldots, N} to 2-d 
//' index \eqn{i = 1, \ldots, a}, \eqn{j = 1, \ldots, n_i}{j = 1, ..., n_i}. 
//' Generally the covariate values from all treatments are stored together in 
//' one vector and \eqn{r = 1, \ldots, N} enumerates the values. For any integer
//' between 1 and \eqn{N}, \code{mapindex()} tells which treatment the 
//' \eqn{r^{th}} value belongs to, and which observation in the identified 
//' treatment.
//' 
//' This is a rewrite of the original \code{mapindex()}function from the PSSIM
//' library using Rcpp to decrease runtime. If you load both this library and 
//' PSSIM, load this one second so the faster functions will overwrite the 
//' slower ones.
//'
//' @param r An integer between 1 and \code{num(n)}.
//' @param n A vector of the sample sizes.
//' 
//' @return A 2-d index, where the first gives which treatment \code{r} belongs 
//' to and the second gives which observation \code{r} is within that treatment.
//' 
//' @examples
//' r = 5; n = c(7, 8); mapindex(r, n)
//' r = 7; n = c(7, 8); mapindex(r, n)
//' r = 9; n = c(7, 8); mapindex(r, n)
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
IntegerVector mapindex(int r, IntegerVector n) {
  Environment base("package:base");
  Function cumsum = base["cumsum"];
  
  IntegerVector sumn = cumsum(n);  
  IntegerVector rem = r - sumn;
  
  int imap = sum(rem > 0); 
  int jmap;
  if (imap < 1) {
    jmap = r - 1;
  } else {
    jmap = rem[imap-1] - 1;
  } 
  // Remove + 1s if you want zero indexing
  return IntegerVector::create(imap + 1, jmap + 1);
}