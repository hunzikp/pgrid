#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getPrioMatrix(IntegerVector gid, NumericVector value) {

  NumericMatrix out(360, 720);
  out.fill(NA_REAL);
  for (int i = 0; i < gid.length(); ++i) {

    int thisGid = gid.at(i) - 1;
    int thisRow = 360 - std::floor(thisGid / 720);
    int thisCol = thisGid - 720*(360 - thisRow);

    out(thisRow, thisCol) = value.at(i);
  }

  return out;
}
