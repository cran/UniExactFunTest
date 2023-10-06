// Cpp code to export
#include <Rcpp.h>
#include "UEFT.h"
#include <string>
using namespace Rcpp;

// format change
vector<vector<int> > matrixtovector(NumericMatrix  input){
  unsigned ncol, nrow;
  ncol = input.ncol();
  nrow = input.nrow();

  vector<vector<int> > ivalue(nrow, vector<int>(ncol));

  for( unsigned i = 0; i < nrow; i++){
    for( unsigned j = 0; j < ncol; j++){
      ivalue[i][j] = input(i,j);
    }
  }
  return ivalue;
}

// Get statistic
// [[Rcpp::export]]
double getFunChisqStat(NumericMatrix  input){
  vector<vector<int> > M = matrixtovector(input);
  return getchisqInteger(M);
}

/////////////////
// UEFT section

// The UEFT code
// [[Rcpp::export]]
double UniEFT(NumericMatrix  input){
  vector<vector<int> > M = matrixtovector(input);
  return  getUniEFT(M, true, false);
}

// The UEFTC code
// [[Rcpp::export]]
double UniEFTC(NumericMatrix  input){
  vector<vector<int> > M = matrixtovector(input);
  return  getUniEFTC(M, false);
}
