# utilities

vecQuadIdx <- function(d){
  # genenrate combinations of seq(d) taken 2 at a time
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
#' @param local.reach approximate size of local neighborhood
#' (treated as Euclidean space)
#'
#' @return a matrix for indices of edge, one row for one edge.
#' @details will based on coordinates, i.e., no edge if coordinates are too far
#' apart, in the sense that difference in coord is greater than \code{range}.
#' @export
#'
#' @examples TBD
allEdge <- function(coord, local.reach){

  stopifnot(all(is.finite(coord)) & all(is.finite(local.reach)))

  if(!is.matrix(coord))
    stop('need more than 1 point as multi-row matrix, check input type.')
  stopifnot(all(local.reach >= 0))
  stopifnot(length(local.reach) == 1 | length(local.reach) == ncol(coord))
  if(length(local.reach) == 1)
    local.reach <- rep(local.reach, ncol(coord))

  res <- locWindow_cpp(
    coord, local.reach, apply(coord, 2, order), apply(coord, 2, rank)
  ) # now it is a list

  # make into a matrix
  tm <- unlist(res)
  mat <- matrix(0, nrow = length(tm), ncol = 2)
  mat[, 1] <- rep(seq_along(res), times = sapply(res, length))
  mat[, 2] <- tm

  return(mat)

#   # index of all possible pairs within range of local.reach
#   mat.order <- apply(coord, 2, order)
#   mat.rank <- apply(coord, 2, rank)
#
#   res <- lapply(seq(nrow(mat.order)), function(i){
#     idx.nbhd <- mat.rank[i, ]
#     idx.nbhd <- pmin(n.pnts, pmax(1, idx.nbhd))
#
#   })
#   idx.pair <- t(utils::combn(nrow(coord.pnts), 2))
#   idx.dist <- idx.dist[
#     sample(nrow(idx.dist), size = 25000),
#   ]

}
