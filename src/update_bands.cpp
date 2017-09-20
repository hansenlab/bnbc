#include <Rcpp.h>
#include <math.h>
#include <cmath>

using namespace Rcpp;


// [[Rcpp::export()]]
NumericMatrix replace(NumericMatrix mat, IntegerMatrix idx, NumericVector update){
  int a,b = 0;
  for (int ii=0; ii < idx.nrow(); ii++){
    a = idx(ii, 0);
    b = idx(ii, 1);
    mat(a, b) = update[ii];
  }
  return mat;
}

// [[Rcpp::export()]]
IntegerMatrix getBandIdxC(int n, int band_no){
  int n_elems = n - band_no + 1;
  IntegerMatrix idx_mat = IntegerMatrix(n_elems, 2);
  int col = band_no-1;
  for (int ii=0; ii < n_elems; ii++){
    idx_mat(ii, 0) = ii;
    idx_mat(ii, 1) = col;
    col++;
  }
  return idx_mat;
}

// [[Rcpp::export()]]
IntegerVector getSeq(int a, int b){
  IntegerVector out = IntegerVector(b - a + 1);
  for (int ii = a; ii <= b; ii++){
    out[ii - a] = ii;
  }
  return out;
}

List deepCopy(List in_list){
  List out_list(in_list.size());
  for (int ii=0; ii < in_list.size(); ii++){
    const NumericMatrix tmp = in_list[ii];
    NumericMatrix tmp2 = Rcpp::clone(tmp);
    out_list[ii] = tmp2;
  }
  return out_list;
}

// [[Rcpp::export()]]
List updateBand(List tact_list, IntegerMatrix idx, NumericMatrix band){
  const List original(tact_list);
  List tact_list2 = deepCopy(original);
  for (int ii=0; ii < band.ncol(); ii++){
    NumericVector band_col = band(_, ii);
    NumericVector v = band_col;
    tact_list2[ii] = replace(tact_list2[ii], idx, band_col);
    NumericMatrix m = tact_list2[ii];
  }
  return tact_list2;
}


// [[Rcpp::export()]]
List updateBands(List tact_list, int n, int bstart, int nbands, List band_mats){
  IntegerVector bands = getSeq(bstart, nbands);
  const List original(tact_list);
  List new_list = deepCopy(original);  // http://gallery.rcpp.org/articles/sorting/ clone
  for (int ii=0; ii < band_mats.size(); ii++){
    IntegerMatrix idx = getBandIdxC(n, bands[ii]);
    new_list = updateBand(new_list, idx, band_mats[ii]);
  }
  return new_list;
}
