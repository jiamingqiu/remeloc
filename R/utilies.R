# utilities

vecQuadIdx <- function(d){
  # generate combinations of seq(d) taken 2 at a time
  # with replacement.
  # Returns: a 2 X d(d-1)/2 where each column is a
  # combination.
  stopifnot(d >= 2)
  return(cbind(
    # diag elements
    t(matrix(seq(d), nrow = d, ncol = 2)),
    utils::combn(d, 2) # non-diag elements
  ))
}

#' Vectorization of terms in quadratic form.
#' @md
#' @param x an array or a matrix
#'
#' @return a matrix with each row being vectorized quadratic form.
#'
#' @details Columns will be in the order of
#'   \eqn{x_1^2, \dots, x_d^2, x_1x_2, \dots, x_{d-1}x_d}.
#' @export
#'
#' @examples TBD
vecQuad <- function(x){
  # vectorization of quadratic form
  # x: an array or a matrix
  # return: a matrix
  if(!is.matrix(x))
    x <- matrix(x, nrow = 1)
  d <- ncol(x)
  if(is.null(colnames(x))){
    colnames(x) <- sprintf('x%s', seq(d))
  }
  # index and names of vectorized result
  idx.col <- vecQuadIdx(d)
  nm.col <- colnames(x)
  nm.col <- sprintf(
    '%s*%s',
    nm.col[idx.col[1, ]],
    nm.col[idx.col[2, ]]
  )
  res <- x[, idx.col[1, ], drop=FALSE] * x[, idx.col[2, ], drop=FALSE]
  colnames(res) <- nm.col
  return(res)
}

#' Quadratic form
#'
#' @param mat.x a matrix
#' @param mat.a symmetric matrix of coefficients.
#'
#' @return a ncol(mat.x) X ncol(mat.x) matrix.
#' @details compute \code{t(mat.x) %*% mat.a %*% mat.x}, so the (i, j) element
#' is \code{mat.x[, i] %*% mat.a %*% mat.x[, j]}.
#' @export
#'
#' @examples TBD
#' @md
quadForm <- function(mat.x, mat.a){
  # compute t(mat.x) %*% mat.a %*% mat.x
  res <- crossprod(mat.x, mat.a)
  res <- tcrossprod(res, t(mat.x))
  return(res)
}

# sigmoid function and its inverse
sigmoid.f <- function(x){
  1 / (1 + exp(-x))
}
sigmoid.f.inv <- function(x){
  log(x / (1 - x))
}

#' Random generator for a local graph
#'
#' @param n number of edge to generate
#' @param coord matrix of coordinates, one row for one point
#' @param local.reach approximate size of local neighborhood
#' (treated as Euclidean space)
#'
#' @return a n-by-2 matrix for indices of edge, one row for one edge.
#' @details will based on coordinates, i.e., no edge if coordinates are too far
#' apart, in the sense that difference in coord is greater than \code{range}.
#' @export
#'
#' @examples TBD
genGraph <- function(n, coord, local.reach){

  if(!is.matrix(coord))
    stop('need more than 1 points.')
  stopifnot(all(local.reach >= 0))
  stopifnot(length(local.reach) == 1 | length(local.reach) == ncol(coord))
  stopifnot(n > 0)

  coord.pnts <- coord

  # index of all possible pairs within range of local.reach

  idx.pair <- t(utils::combn(nrow(coord.pnts), 2))
  idx.dist <- idx.dist[
    sample(nrow(idx.dist), size = 25000),
  ]

  return(n * range)
}


#' get all possible local edge
#'
#' @param coord matrix of coordinates, one row for one point
#' @param local.reach approximate size of local neighborhood square
#' (treated as Euclidean space), recycled if necessary
#'
#' @return a 2-cols matrix for indices of edge, one row for one edge.
#' @details based on coordinates, no edge if coordinates are too far apart, in
#' the sense that difference in coord is greater than \code{local.reach}.
#' @export
#'
#' @examples
#' set.seed(1)
#' d <- 3
#' coord <- matrix(runif(d * 50), ncol = d)
#' edge <- allEdge(coord, 0.25)
#' plotGraph(coord, edge)
#' # another example
#' coord <- matrix(0, nrow = 5, ncol = 2)
#' coord[-1, ] <- as.matrix(expand.grid(v1 = c(-1, 1), v2 = c(-1, 1)))
#' edge <- allEdge(coord, 1)
#' plotGraph(coord, edge)
allEdge <- function(coord, local.reach){
  # a function finding neighbours
  # args:
  #   coords: a matrix for coordinates of points, 1 row = 1 point
  #   width: width of intervals
  #   return: list(default)/data.frame, what to return
  # returns: a list/data.frame of 2 column with elements being seq(nrow(coords)),
  #   for example, a row of (1, 2) means the point coord[2, ] is within
  #   the box of width centering at point coord[1, ].
  # Details: only judging by difference in coordinates, does not care about
  #   geometry. Also, the ith list will start from i.

  # coord <- matrix(rnorm(2 * 5000), ncol = 2, nrow = 5000)
  # local.reach <- 0.15
  stopifnot(all(is.finite(coord)) & all(is.finite(local.reach)))

  if(!is.matrix(coord))
    stop('need more than 1 point as multi-row matrix, check input type.')
  stopifnot(all(local.reach >= 0))
  stopifnot(length(local.reach) == 1 | length(local.reach) == ncol(coord))
  if(length(local.reach) == 1)
    local.reach <- rep(local.reach, ncol(coord))

  d <- ncol(coord)
  n <- nrow(coord)
  stopifnot(n > 1)

  if(n > 1000){
    # larger n, R is faster
    tm.res <- lapply(seq(n - 1), function(idx){
      c(idx, which(
        colSums(
          abs(t(coord[seq(idx + 1, n), ]) - coord[idx, ]) <= local.reach
        ) == d
      ) + idx)
    })
    tm.res[n] <- n # for completeness
  }else{
    # smaller, rcpp faster, very weird...
    tm.res <- locWindow_cpp(
      coord, local.reach,
      apply(coord, 2, order), apply(coord, 2, rank, ties.method = 'first')
    ) # now it is a list
  }

  # make into a matrix
  tm <- unlist(tm.res)
  mat <- matrix(0L, nrow = length(tm), ncol = 2) # 0L to keep integer type
  mat[, 1] <- rep(seq_along(tm.res), times = sapply(tm.res, length))
  mat[, 2] <- tm

  return(mat)
}

# allEdgeCpp <- function(coord, local.reach){
#
#   stopifnot(all(is.finite(coord)) & all(is.finite(local.reach)))
#
#   if(!is.matrix(coord))
#     stop('need more than 1 point as multi-row matrix, check input type.')
#   stopifnot(all(local.reach >= 0))
#   stopifnot(length(local.reach) == 1 | length(local.reach) == ncol(coord))
#   if(length(local.reach) == 1)
#     local.reach <- rep(local.reach, ncol(coord))
#
#   res <- locWindow_cpp(
#     coord, local.reach,
#     apply(coord, 2, order), apply(coord, 2, rank, ties.method = 'first')
#   ) # now it is a list
#
#   # browser();QWER
#   # make into a matrix
#   tm <- unlist(res)
#   mat <- matrix(0L, nrow = length(tm), ncol = 2) # 0L to keep integer type
#   mat[, 1] <- rep(seq_along(res), times = sapply(res, length))
#   mat[, 2] <- tm
#
#   # # drop edge pointing to self
#   # mat <- mat[mat[, 1] != mat[, 2], ]
#
#   return(mat)
#
#   #   # index of all possible pairs within range of local.reach
#   #   mat.order <- apply(coord, 2, order)
#   #   mat.rank <- apply(coord, 2, rank)
#   #
#   #   res <- lapply(seq(nrow(mat.order)), function(i){
#   #     idx.nbhd <- mat.rank[i, ]
#   #     idx.nbhd <- pmin(n.pnts, pmax(1, idx.nbhd))
#   #
#   #   })
#   #   idx.pair <- t(utils::combn(nrow(coord.pnts), 2))
#   #   idx.dist <- idx.dist[
#   #     sample(nrow(idx.dist), size = 25000),
#   #   ]
#
# }
