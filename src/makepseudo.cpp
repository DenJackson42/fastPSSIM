#include "RcppArmadillo.h"
using namespace Rcpp;

// Helper function for makepseudo, don't export to R
List create_newtrt(int i, int j, IntegerVector n, NumericMatrix alltrt) {
  
  // Get some R functions that will be useful later
  Environment base("package:base");
  Function cumsum = base["cumsum"];
  Function rank = base["rank"]; 
  
  // Get cumulative sums for each treatment level, to be used for indexing later
  IntegerVector cum_n = cumsum(n);
  
  NumericMatrix newtrt;
  int total = 0, jj = 0;
  
  // If we are in the first treatment...
  if (i == 1) {
    // If j is an obs. within the first treatment...
    if (j <= n[0]) {
      // newtrt is initalized as the entire first treatment
      newtrt = alltrt(_, Range(0, n[0] - 1));
      total = newtrt.ncol();
      jj = j;
    } 
    // If j is not in the first treatment...
    if (j > n[0]) {
      // Get all obs. in first treatment
      newtrt = alltrt(_, Range(0, n[0] - 1));
      // Get column of jth obs.
      NumericVector col = alltrt(_, j - 1);
      newtrt = cbind(newtrt, col);
      total = newtrt.ncol();
      jj = total;
    } 
  } 
  if (i > 1) { // If we aren't in the first treatment...
    // Get upper and lower bounds of our treatment
    int trt_upper = cum_n[i - 1];// upper 
    int trt_lower = cum_n[i - 2]; // lower
    
    if (j > trt_lower && j <= trt_upper) { // If jth obs. is in ith trt...
      newtrt = alltrt(_, Range(trt_lower, trt_upper - 1));
      total = newtrt.ncol();
      jj = j - trt_lower; // Gets index of jth obs within its treatment
    } else { // If jth obs. is not in the ith treatment...
      // Get all obs. in the ith treatment
      newtrt = alltrt(_, Range(trt_lower, trt_upper - 1));
      // Get column of jth obs.
      NumericVector col = alltrt(_, j - 1);
      newtrt = cbind(newtrt, col);
      total = newtrt.ncol();
      jj = total;
    }
  }
  // Make sure ranks of X are correct, which are in 3rd row of newtrt
  NumericVector ranks = rank(newtrt.row(1));
  newtrt(2, _) = ranks;
  
  // Create flag variable
  LogicalVector flag(total);
  for (int idx = 0; idx < total; idx++) {
    flag[idx] = (jj == total && jj > n[i - 1] && idx != jj - 1) || 
      (jj <= n[i - 1]);
  }
  
  // Get the jjth obs' rank
  int target_rank = newtrt(2, jj - 1);
  
  // If jj is the total # of obs in the trt and it is past the number of obs 
  // in the ith trt...
  if (jj == total && jj > n[i - 1]) {
    // Remove the jjth column and rerank the obs, and bring total down 1
    NumericVector temp = newtrt.row(1);
    temp.erase(temp.begin() + jj - 1);
    NumericVector temp_ranks = rank(temp);
    newtrt.row(2) = temp_ranks;
    total--;
  }
  
  // Filter new_trt to only grab columns where flag = TRUE
  NumericMatrix final_newtrt(newtrt.nrow(), sum(flag));
  int col = 0;
  for (int idx = 0; idx < total; idx++) {
    if (flag[idx]) {
      final_newtrt(_, col) = newtrt(_, idx);
      col++;
    }
  }
  
  return List::create(Named("newtrt") = final_newtrt,
                      Named("target_rank") = target_rank);
}

//Helper function for makepseudo, don't export to R
IntegerVector find_nearest_neighbors(int target_rank, int k,
                                     NumericMatrix newtrt) {
  
  // Get the order function from R. Might code own, but it's slower for large ns
  Environment base("package:base");
  Function order = base["order"];
  
  IntegerVector ordk(k);
  int total = newtrt.ncol(); 
  
  if (target_rank <= ((k - 1)/2)) {
    ordk = order(newtrt(2, _));
    ordk = ordk[Range(0, k-1)]; 
  }
  if (target_rank > (total - ((k-1)/2))) {
    ordk = order(total - newtrt(2, _));
    ordk = ordk[Range(0, k - 1)];
  }
  if (target_rank > (k-1)/2 && target_rank <= (total - (k-1)/2)) {
    ordk = order((abs(newtrt(2, _) - target_rank)));
    ordk = ordk[Range(0, k - 1)];
  }
  return ordk;
}

//' Nearest neighbor augmentation based on ranks
//' 
//' @description
//' The function \code{makepseudo()} performs the nearest neighbor augmentation 
//' based on the rank of covariate values according to the scheme described on 
//' page 410-411 of Wang, Tolos and Wang (2010). 
//' 
//' This is a rewrite of the original \code{makepseudo()} function from the 
//' PSSIM library using Rcpp to decrease runtime. If you load both this library 
//' and PSSIM, load this one second so the faster functions will overwrite the 
//' slower ones.
//'
//' @param N Total number of covariate values.
//' @param n Vector of sample sizes from all treatments.
//' @param k Number of nearest neighbors.
//' @param a Number of treatment levels in the data.
//' @param alltrt A matrix of dimension 3x\code{N}, with first two rows being Y 
//'   and X, and a third row with the rank of the X values within the same 
//'   treatment level. 
//' 
//' @return A list containing the following:
//' \itemize{
//'   \item \strong{psudo}: A 3D array of dimensions (\code{k}, \code{a}, 
//'   \code{N}) that stores the augmented observations based on the k-nearest 
//'   neighbor rule in Wang, Tolos and Wang (2010).
//'   \item \strong{index}: A 3D array of dimensions (\code{k}, \code{a}, 
//'   \code{N}) that stores the indices of the observations used for 
//'   augmentation.
//' }
//' 
//' @examples
//'  a = 2; n = c(7, 9); N = sum(n);  X = runif(N);
//'  trt = c(rep(1, n[1]), rep(2, n[2])); e = rnorm(N, 0, 0.1)
//'  Y = ifelse(trt == 1, 4*(X-0.5)^2+e, 2*X+e)
//'  ranksuse = unlist(tapply(X, trt, rank))
//'  alltrt = rbind(Y, X, ranksuse)
//'  aug = makepseudo(N, n, k=3, a, alltrt)
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
List makepseudo(int N, IntegerVector n, int k, int a, NumericMatrix alltrt) {
  NumericVector psudo(a * N * k);
  NumericVector index(a * N * k);
  
  // This puts the 1d matrices into 3d matrices when returned  
  IntegerVector dim = IntegerVector::create(a, N, k);
  psudo.attr("dim") = dim;
  index.attr("dim") = dim;
  
  // Formula to access [i, j, k]th element:
  // index = (i-1) + (j-1)*a + (k-1)*a*N
  // The - 1s are because C++ uses 0 indexing
  
  // Fills psudo and index with pseudo-replicates
  for (int i = 1; i <= a; i++) {
    for (int j = 1; j <= N; j++) {
      // First we get the matrix newtrt, which returns two objects:
      // newtrt, a NumericMatrix with the Y1, X, and Rank of all the obs. in
      // the ith treatment
      // target_rank, the rank of the jth obs in its treatment group (type: int)
      List newtrt_res = create_newtrt(i, j, n, alltrt);
      int target_rank = as<int>(newtrt_res["target_rank"]);
      NumericMatrix newtrt = as<NumericMatrix>(newtrt_res["newtrt"]);
      
      // Now get the indexes of the jth observation's pseudo-replicates
      IntegerVector ordk = find_nearest_neighbors(target_rank, k, newtrt);
      
      // Fill psudo and index:
      // Get the starting index
      int stIdx = (i-1) + (j-1) * a;
      
      for (int l = 0; l < k; l++) {
        // Get index for the [i, j, (k-1)]th element
        int idx = stIdx + l*a*N;
        psudo[idx] = newtrt(0, ordk[l] - 1); // - 1 because C++ has 0 indexing
        index[idx] = ordk[l]; // Add a - 1 for 0 indexing if so desired
      }
    }
  }
  
  return List::create(Named("psudo") = psudo,
                      Named("index") = index);
}