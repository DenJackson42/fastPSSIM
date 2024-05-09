#include "RcppArmadillo.h"
using namespace Rcpp;

List makepseudo(int N, IntegerVector n, int k, int a, NumericMatrix alltrt);
IntegerVector trt_position(int i1, IntegerVector n);
IntegerVector mapindex(int r, IntegerVector n);

//' Test of independence in the presence of heteroscedastic treatments
//' 
//' @description
//' \code{NPtest_indept()} performs the test of independence between the 
//' response variable and a single covariate when there is potentially 
//' heteroscedastic treatment effects present (see Wang, Tolos and Wang (2010)).
//' 
//' This is a rewrite of the original \code{NPtest_indept()} function from the
//' PSSIM library using Rcpp to decrease runtime. If you load both this library 
//' and PSSIM, load this one second so the faster functions will overwrite the 
//' slower ones.
//'
//' @param dat A dataframe with three columns named X, trt, and Y, where X is
//' the covariate, trt is the treatment level, and Y is the response variable.
//' @param k An odd integer to specify the number of nearest neighbors to be
//' used in augmentation. Generally recommend to use 3, 5, or 7.
//' 
//' @return A list containing the following variables:
//' \itemize{
//'   \item \strong{Asys_var}: The asymptotic variance for the test statistics.
//'   \item \strong{Tstat}: The test statistic.
//'   \item \strong{pvalue}: The p-value of the test under \eqn{H_0}:
//'   There is independence between X and Y.
//' }
//' 
//' @examples
//' n = 64;  X = runif(n); trt = gl(2, n/2)
//' e = rnorm(n, 0, 0.1)
//' Y = ifelse(trt == 1, 4*(X-0.5)^2+e, 2*X+e)
//' dat = data.frame(X, Y, trt)
//' NPtest_indept(dat, k = 7)
//'  
//' @seealso \code{\link[=makepseudo]{makepseudo}}, 
//' \code{\link[=mapindex]{mapindex}}, and  
//' \code{\link[=trt_position]{trt_position}} are all called within this.
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
List NPtest_indept(DataFrame dat, int k = 7) {
  
  Environment base("package:base");
  Function rank = base["rank"]; 
  Function order = base["order"]; 
  
  NumericVector X = dat["X"];
  IntegerVector trt = dat["trt"];
  NumericVector Y1 = dat["Y"];
  
  double var_Y1 = var(Y1);
  int N = X.size();
  
  if (var_Y1 > 1e-6) {
    IntegerVector trt_lvls = unique(trt).sort();
    int a = trt_lvls.length();
    IntegerVector ranksuse(N);
    
    // To keep track of number of obs in each trt
    IntegerVector n(2);
    
    // Loop through each treatment and rank X values, this fills ranksuse
    for (int i = 0; i < a; i++) {
      //int lvl = trt_lvls[i];
      // Get index for Xs in current trt
      LogicalVector idx = (trt == trt_lvls[i]);
      NumericVector subX = X[idx];
      IntegerVector ranked = rank(subX);
      
      // Place ranked values back into the appropriate positions in ranksuse
      int k = 0;
      for (int j = 0; j < N; j++) {
        if (idx[j]) {
          ranksuse[j] = ranked[k];
          k++;
        }
      }
      // At this point, k should be # of obs in the current trt, so add to n
      n[i] = k;
      
    }
    // Now rbind Y1, X, and ranksuse together
    NumericMatrix alltrt(3, N);
    alltrt(0, _) = Y1;
    alltrt(1, _) = X;
    alltrt(2, _) = ranksuse;
    
    // This rearranges the data in alltrt
    // a = # of treatments
    for (int i1 = 1; i1 <= a; i1++) {
      IntegerVector locationi1 = trt_position(i1, n);
      int lower = locationi1[0] - 1;
      int upper = locationi1[1] - 1;
      // Have to extract row first and then subset it
      NumericVector tempRow = alltrt(1, _);
      IntegerVector orderwant = order(tempRow[Range(lower, upper)]);
      orderwant = orderwant + lower - 1;
      // Now reorder the columns in alltrt
      NumericMatrix temp = clone(alltrt); // Clone alltrt to use as a ref.
      for (int col = lower; col <= upper; col++) {
        alltrt(_, col) = temp(_, orderwant[col - lower]);
      }
    }
    
    List psudodat = makepseudo(N, n, k, a, alltrt);
    NumericVector psudo = psudodat["psudo"]; 
    NumericVector index = psudodat["index"];
    
    // Now we need to get the cellmean, colmean, and sig from psudo
    // Initialize holders:
    arma::mat cellmean(a, N);
    NumericVector colmean(N);
    // For matrix with \hat\sigma_i^2(X_{ij}) = sigXij[i, j]
    NumericMatrix sigXij(a, N);
    // To hold mean of each row in psudo
    NumericVector meanrk(a);
    
    // Calculate cellmean (takes mean of all row 1 column 1 values for all k
    // matrices, then all row 1 column 2 values, and so on)
    for (int i = 0; i < a; i++) {
      // Holds amount to get rowSum later
      double rowSum = 0;
      for (int j = 0; j < N; j++) {
        
        double sum = 0;
        // Hold values to take variance of later
        NumericVector tempVals(k);
        
        for (int l = 0; l < k; l++) {
          // i + j * nrow + l*nrow*ncol for indexing 3d matrix
          double val = psudo[i + j * a + l * a * N];
          sum += val;
          tempVals[l] = val;
          rowSum += val;
        }
        cellmean(i, j) = sum / k;
        sigXij(i, j) = var(tempVals);
      }
      meanrk[i] = rowSum / (N*k);
    }
    
    // Calculate colmean (mean of all values in column 1, 2, etc. across all k)
    for (int j = 0; j < N; j++) {
      double sum = 0;
      for (int i = 0; i < a; i++) {
        // Averaging the cellmean columns gets us the column averages
        sum += cellmean(i, j);
      }
      colmean(j) = sum / a;
    }
    // Using arma::mat because it's just a lot easier for covariance
    arma::mat cellmeanT = arma::trans(cellmean);
    arma::mat sig = arma::cov(cellmeanT);
    // diagonal part gives the \hat\sigma_{1,i}^2
    // and off-diagonal part gives \hat\sigma_{1,i_1, i_2}
    
    // Now we need to calculate MSTphi
    double sumOfSquares = 0.0;
    
    // Get sum of the squared differences
    for (int i = 0; i < a; i++) {
      for (int j = 0; j < N; j++) {
        double diff = cellmean(i, j) - meanrk[i];
        sumOfSquares += diff * diff;
      }
    }
    // Calculate MSTphi
    double MSTphi = k * sumOfSquares / ((N - 1) * a);
    
    // Now to get MSE
    sumOfSquares = 0.0; // Reset sumOfSquares
    
    // Take the difference between each value in psudo and its corresponding
    // cellmean value
    // if (use_parallel) {
    //   #pragma omp parallel for reduction(+:sumOfSquares)
    //   for (int i = 0; i < a; i++) {
    //     for (int j = 0; j < N; j++) {
    //       for (int l = 0; l < k; l++) {
    //         // To index a 3d matrix, idx = i + j * nrow + l * nrow * ncol
    //         double diff = psudo[i + j * a + l * a * N] - cellmean(i, j);
    //         sumOfSquares += diff * diff;
    //       }
    //     }
    //   }
    // } else {}
    for (int i = 0; i < a; i++) {
      for (int j = 0; j < N; j++) {
        for (int l = 0; l < k; l++) {
          // To index a 3d matrix, idx = i + j * nrow + l * nrow * ncol
          double diff = psudo[i + j * a + l * a * N] - cellmean(i, j);
          sumOfSquares += diff * diff;
        }
      }
    }
    // Calculate MSE
    double MSE = sumOfSquares / (N*a*(k-1));
    
    double Tss = (sqrt(N) * (MSTphi - MSE));
    
    // Calculate estimate of variance for test statistics
    // count is a matrix; first three columns give the value of i1, j2, i;
    //       the last column gives the number of times X_{ij_2} is used in
    //       construction of windows for all covariate values in group i_1
    
    // Initialize matrix to hold results
    NumericMatrix count(a*a*N, 4);
    
    // Convert index to a cube because I'm not smart enough to get it working as
    // as flat NumericVector
    arma::cube indexCube(index.begin(), a, N, k);
    
    int whereini = -1; // Start from -1 since Rcpp is zero-index
    
    for (int i1 = 1; i1 <= a; i1++) {
      IntegerVector locationi1 = trt_position(i1, n);
      int lower = locationi1[0] - 1;
      int upper = locationi1[1] - 1;
      // Don't think parallel will work because of whereini
      for (int j2 = 1; j2 <= N; j2++) {
        for (int i = 0; i < a; i++) {
          whereini++;
          IntegerVector whereisXij2 = mapindex(j2, n);
          
          int counti1j2i = 0;
          int targetRank = whereisXij2[1] - 1;
          // Had to separate the booleans because 0 is an actual index in Rcpp
          bool isCorrectTrt = (whereisXij2[0] - 1) == i;
          
          for (int j = lower; j <= upper; j++) {
            for (int l = 0; l < k; l++) {
              if (isCorrectTrt && (indexCube(i, j, l) - 1) == targetRank) {
                counti1j2i++;
              }
            }
          }
          
          NumericMatrix::Row row = count.row(whereini);
          row[0] = i1-1;
          row[1] = j2;
          row[2] = i;
          row[3] = counti1j2i;
        }
      }
    }
    
    // Get subcount, the entries in the count where i1 != i3
    // To store rows temporarily, we don't know how many so this is dynamic
    std::vector<NumericVector> temp;  
    
    for (int i = 0; i < count.nrow(); i++) {
      if (count(i, 0) != count(i, 2)) {
        temp.push_back(count.row(i));  // Store row in temp
      }
    }
    
    // Now create subcount from temp
    int tempL = temp.size();
    NumericMatrix subcount(tempL, count.ncol());
    for (int i = 0; i < tempL; i++) {
      subcount.row(i) = temp[i];
    }
    
    // Using a map to emulate the behavior of tapply with sum (thanks ChatGPT)
    std::map<std::pair<int, int>, double> sums;
    
    for (int i = 0; i < subcount.nrow(); i++) {
      int j2 = subcount(i, 1);  
      int i_ = subcount(i, 2);
      double value = subcount(i, 3);
      
      // Create a pair for the map key
      std::pair<int, int> key = std::make_pair(j2, i_);
      
      // Add to the sum for this key
      sums[key] += value;
    }
    
    // Prepare output equivalent to tapply output, excpet as a vector not a matrix
    NumericVector prodcount1_1(sums.size());
    int idx = 0;
    for (auto& kv : sums) {
      prodcount1_1[idx++] = kv.second;
    }
    
    double tau3 = 0;
    
    // Trying to parallelize this loop increases runtime across the board
    //#pragma omp parallel for reduction(+:tau3)
    for (int i = 0; i < a; i++) {
      // Get start position of the current trt
      int starti = trt_position(i+1, n)[0] - 1;
      for (int jp = starti + 1; jp < starti+n[i]; jp++) { // Skip first obs
        for (int j = std::max(0, (jp-k+1)); j < jp; j++) {
          // Math I don't get :( p.414 of Wang
          // j * ncol + i to index flattened 2d array
          double Bijjp = (prodcount1_1[j * a + i]/k + 1)*(k-jp+j)*(jp-j<=k-1);
          double tauAdd = ((Bijjp*Bijjp + Bijjp) - 2*(jp-j<=((k-1)/2)))*
            (jp-j<=k-1) * sigXij(i, j) * sigXij(i, jp) * (j!=jp);
          tau3 += tauAdd;
        }
      }
    }
    // This completes the formula for (gamma-hat_)N^2 on p. 414, middle of page
    tau3 = tau3*4/(N*a*a*(k-1)*(k-1));
    
    double tauAsys = tau3;
    
    // Get p-value using test statistic Tss / sqrt(tauAsys)
    double pvalue_sim = R::pnorm(Tss/sqrt(tauAsys), 0, 1, 0, 0);
    
    return List::create(Named("Asys_var") = tauAsys,
                        Named("Tstat") = Tss,
                        Named("pvalue") = pvalue_sim);
    
  } else {
    // If the var is less than the number basically just give up, do not reject H0
    int Tss = 0; int tauAsys = 1e+8;  int pvalue_sim = 1;
    return List::create(Named("Asys_var") = tauAsys,
                        Named("Tstat") = Tss,
                        Named("pvalue") = pvalue_sim);
  }
}