// [[Rcpp::plugins("cpp17")]]
#include <cstdio>

// edlib functions dependency
#include "edlib/edlib.h"
#include <Rcpp.h>

using namespace Rcpp;

// fast levenshtein distance calculation with edlib 
int Levenshtein(const char* seq_1, const char* seq_2){
  // define input sequences as string constants
  static std::string work1; static std::string work2;
  work1 = seq_1; work2 = seq_2;
  
  const char *query  = work1.c_str();
  const char *target = work2.c_str();
  
  // edlib function
  EdlibAlignResult result = edlibAlign(
    query, strlen(query), target, strlen(target), 
    edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0)
  );
  
  // return edit distance
  if(result.status == EDLIB_STATUS_OK){
    return result.editDistance;
  }

  return 0;
}

// [[Rcpp::export]]
NumericMatrix LevenshteinLoop(Rcpp::CharacterMatrix mat){
  int ncol = mat.ncol()-1;
  int nrow = mat.nrow();

  // initialise numeric matrix for storing results
  Rcpp::NumericMatrix out(nrow, ncol);

  for(int i = 0; i < nrow; ++i){
    for(int j = 0; j < ncol; ++j){
      out(i, j) = Levenshtein(mat(i, j), mat(i, ncol));
    }
  }

  return out;
}