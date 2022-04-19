#include <Rcpp.h>
using namespace Rcpp;

// Functions used to determine local graph
// essentially index all possible local edges

// function declaration
IntegerVector orderCpp(NumericVector x);
IntegerVector orderStable(NumericVector x);
IntegerVector rankCpp(NumericVector x);

// [[Rcpp::export]]
List locWindow_cpp(
    NumericMatrix coord, NumericVector local_reach,
    IntegerMatrix marginOrder, IntegerMatrix marginRank
) {

  // take matrix of coordinates (coord, 1row = 1point)
  // return index of all possible local edges, where "local" defined
  // as difference of coord <= local_reach.
  // the marginOrder is the order of each column of coord, 1 to n
  // the marginRank is the rank of each column of coord, 1 to n
  // the marginOrder and marginRank should start with 1.

  int n = coord.nrow();
  int p = coord.ncol();

  if (
      Rcpp::min(marginOrder) != 1 || Rcpp::max(marginOrder) != n ||
        marginOrder.nrow() != n || marginOrder.ncol() != p
    )
    Rcpp::stop("check input marginOrder");
  if (
      Rcpp::min(marginRank) != 1 || Rcpp::max(marginRank) != n ||
        marginRank.nrow() != n || marginRank.ncol() != p
  )
    Rcpp::stop("check input marginRank");
  if (Rcpp::min(local_reach) < 0 || local_reach.size() != p)
    Rcpp::stop("check input local_reach");

  marginOrder = marginOrder - 1;
  marginRank = marginRank - 1;

  // Nullable<IntegerMatrix> marginOrder = R_NilValue,
  // Nullable<IntegerMatrix> marginRank = R_NilValue
  // if (local_reach.size() != p) {
  //   NumericVector use_reach(p);
  //   std::fill(use_reach.begin(), use_reach.end(), local_reach);
  // } else {
  //   NumericVector use_reach = local_reach;
  // }
  // // get sort index if not provided
  // if (marginOrder.isNull()) {
  //   IntegerMatrix marginOrder(n, p);
  //   for (int j = 0; j < p; j++) {
  //     marginOrder(_, j) = orderCpp(coord(_, j));
  //   }
  // }
  // if (marginRank.isNull()) {
  //   IntegerMatrix marginRank(n, p);
  //   for (int j = 0; j < p; j++) {
  //     marginRank(_, j) = rankCpp(coord(_, j));
  //   }
  // }

  // figure out local neighborhood
  List res(n);
  // LogicalMatrix mat_nbhd(n, p);
  // intersect bands
  LogicalVector tm_nbhd(n);
  IntegerVector tm_idx(n); // 1 to n for index
  for (int i = 0; i < n; i++) {
    tm_idx[i] = i + 1;
  }
  bool keep_pnt;

  for (int i = 0; i < n; i++) {

    if (i % 50000 == 0) { // for interrupt in case n very large
      Rcpp::checkUserInterrupt();
    }

    //  for the i-th point, for unique pairs, only consider larger than i
    std::fill(tm_nbhd.begin(), tm_nbhd.begin() + i, FALSE);
    std::fill(tm_nbhd.begin() + i, tm_nbhd.end(), TRUE);

    // get bands per j-th coord and overlay them
    for (int j = 0; j < p; j++) {

      tm_nbhd[i] = TRUE;
      int rankThisPnt = marginRank(i, j);

      // search to the larger side
      bool inrange = TRUE;
      for (int rankR = rankThisPnt; rankR < n; rankR++) {
        int pR = marginOrder(rankR, j); //pointer to the right point
        if (tm_nbhd[pR] && inrange) { // if still in range
          keep_pnt =
            std::abs(coord(i, j) - coord(pR, j)) <= local_reach[j];
          tm_nbhd[pR] = keep_pnt && tm_nbhd[pR];
          if (!keep_pnt) { // flag outside local range
            inrange = FALSE;
          }
        } else { // already outside
          tm_nbhd[pR] = FALSE;
        }
      }
      // search to the smaller side
      inrange = TRUE;
      for (int rankL = rankThisPnt; rankL >= 0; rankL--) {
        int pL = marginOrder(rankL, j); //pointer to the right point
        if (tm_nbhd[pL] && inrange) { // if still in range
          keep_pnt =
            std::abs(coord(i, j) - coord(pL, j)) <= local_reach[j];
          tm_nbhd[pL] = keep_pnt && tm_nbhd[pL];
          if (!keep_pnt) { // flag outside local range
            inrange = FALSE;
          }
        } else { // already outside
          tm_nbhd[pL] = FALSE;
        }
      }

    }

    res[i] = tm_idx[tm_nbhd];

  }

  return res;

}

// no need to export the following, since the R's order & rank is faster now
IntegerVector orderCpp(NumericVector x) {
  // return order of x
  // https://stackoverflow.com/questions/21609934/ordering-permutation-in-rcpp-i-e-baseorder
  // this warning makes things very slow for large array
  // if (is_true(any(duplicated(x)))) {
  //   Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  // }
  NumericVector sorted = clone(x);
  std::sort(sorted.begin(), sorted.end());
  return match(sorted, x);
}

IntegerVector orderStable(NumericVector x) {
  // return order of x
  // https://stackoverflow.com/questions/21609934/ordering-permutation-in-rcpp-i-e-baseorder
  // if (is_true(any(duplicated(x)))) {
  //   Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  // }
  NumericVector sorted = clone(x);
  std::stable_sort(sorted.begin(), sorted.end());
  return match(sorted, x);
}

IntegerVector rankCpp(NumericVector x) {
  // return rank of x
  NumericVector sorted = clone(x);
  std::sort(sorted.begin(), sorted.end());
  return match(x, sorted);
}
