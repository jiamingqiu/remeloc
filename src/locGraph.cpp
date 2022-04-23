#include <Rcpp.h>
using namespace Rcpp;

// Functions used to determine local graph
// essentially index all possible local edges

// function declaration
IntegerVector orderCpp(NumericVector x);
IntegerVector orderStable(NumericVector x);
IntegerVector rankCpp(NumericVector x);

// [[Rcpp::export]]
List allEdge_cpp(
    NumericMatrix coord, NumericVector local_reach,
    IntegerMatrix marginOrder, IntegerMatrix marginRank
) {

  // take matrix of coordinates (coord, 1row = 1point)
  // return index of all possible local edges, where "local" defined
  // as difference of coord <= local_reach.
  // the marginOrder is the order of each column of coord, 1 to n
  // the marginRank is the rank of each column of coord, 1 to n
  // the marginOrder and marginRank should start with 1.

  R_xlen_t n = coord.nrow();
  int d = coord.ncol();

  // for the sake of safety
  if (
      Rcpp::min(marginOrder) != 1 || Rcpp::max(marginOrder) != n ||
        marginOrder.nrow() != n || marginOrder.ncol() != d
    )
    Rcpp::stop("check input marginOrder");
  if (
      Rcpp::min(marginRank) != 1 || Rcpp::max(marginRank) != n ||
        marginRank.nrow() != n || marginRank.ncol() != d
  )
    Rcpp::stop("check input marginRank");
  if (Rcpp::min(local_reach) < 0 || local_reach.size() != d)
    Rcpp::stop("check input local_reach");

  marginOrder = marginOrder - 1;
  marginRank = marginRank - 1;

  // Nullable<IntegerMatrix> marginOrder = R_NilValue,
  // Nullable<IntegerMatrix> marginRank = R_NilValue
  // if (local_reach.size() != d) {
  //   NumericVector use_reach(d);
  //   std::fill(use_reach.begin(), use_reach.end(), local_reach);
  // } else {
  //   NumericVector use_reach = local_reach;
  // }
  // // get sort index if not provided
  // if (marginOrder.isNull()) {
  //   IntegerMatrix marginOrder(n, d);
  //   for (int j = 0; j < d; j++) {
  //     marginOrder(_, j) = orderCpp(coord(_, j));
  //   }
  // }
  // if (marginRank.isNull()) {
  //   IntegerMatrix marginRank(n, d);
  //   for (int j = 0; j < d; j++) {
  //     marginRank(_, j) = rankCpp(coord(_, j));
  //   }
  // }

  // figure out local neighborhood
  List res(n);
  // intersect bands
  LogicalVector tm_nbhd(n);
  LogicalVector loc_nbhd(n);
  IntegerVector tm_idx(n); // 1 to n for index
  for (R_xlen_t i = 0; i < n; i++) {
    tm_idx[i] = i + 1;
  }
  bool keep_pnt;

  for (R_xlen_t i = 0; i < n; i++) {

    if (i % 7500 == 0) { // for interrupt in case n very large
      Rcpp::checkUserInterrupt();
    }

    //  for the i-th point, for unique pairs, only consider larger than i
    std::fill(loc_nbhd.begin(), loc_nbhd.begin() + i, FALSE);
    std::fill(loc_nbhd.begin() + i, loc_nbhd.end(), TRUE);

    // get bands per j-th coord and overlay them
    for (int j = 0; j < d; j++) {

      // vector recording per-j (col)
      std::fill(tm_nbhd.begin(), tm_nbhd.end(), FALSE);

      R_xlen_t rankThisPnt = marginRank(i, j);
      R_xlen_t rankR = rankThisPnt; // rank of current pointer to larger side
      R_xlen_t rankL= rankThisPnt; // rank of current pointer to smaller side

      // search to the larger side
      while (rankR < n) {
        R_xlen_t pR = marginOrder(rankR, j); //pointer to the right point
        keep_pnt =
          std::abs(coord(i, j) - coord(pR, j)) <= local_reach[j];
        // mat_nbhd(pR, j) = keep_pnt;
        tm_nbhd[pR] = keep_pnt;
        if (keep_pnt) {
          rankR = rankR + 1; // still inside the window, next
        } else {
          rankR = n + 1; // out of range already
        }
      }
      // search to the smaller side
      while (rankL >= 0) {
        R_xlen_t pL = marginOrder(rankL, j); //pointer to the right point
        keep_pnt =
          std::abs(coord(i, j) - coord(pL, j)) <= local_reach[j];
        // mat_nbhd(pL, j) = keep_pnt;
        tm_nbhd[pL] = keep_pnt;
        if (keep_pnt) {
          rankL = rankL - 1; // still inside the window, next
        } else {
          rankL = -1; // out of range already
        }
      }

      // overlay (&&) to local
      for (R_xlen_t ii = i; ii < n; ii++) {
        loc_nbhd[ii] = loc_nbhd[ii] && tm_nbhd[ii];
      }
      // loc_nbhd[i] = TRUE;

    }

    res[i] = tm_idx[loc_nbhd];

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
