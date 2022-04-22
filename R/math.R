### math functions


##### model spaces: Euclidean, Sphere, and Poincare disk #######################

# what is needed:
#   coordinate system and random generator
#   corresponding metric tensor and geodesic distance
#   embedding into Euclidean space
#   (maybe) exponential/logarithmic map
spaceEuclidean <- function(d){
  # d-dimensional Euclidean space under Cartesian coordinates

  stopifnot(d >= 1)

  genPnt <- function(n){
    # random point generator, returns a matrix, one row = one point
    stopifnot(n > 0)
    return(matrix(runif(n * d), ncol = d, nrow = n))
  }

  # metric tensor
  metric <- function(x = NULL){
    return(diag(rep(1, d)))
  }
  # geodesic distance
  dist <- function(x, y){
    # function to compute geodesic distance on Euclidean space
    # x, y in cartesian coord, one row for one point

    if(!is.matrix(x)) x <- matrix(x, nrow = 1)
    if(!is.matrix(y)) y <- matrix(y, nrow = 1)
    stopifnot(all(dim(x) == dim(y)))

    return(
      sqrt(rowSums((x - y) ^ 2))
    )
  }

  return(list(
    genPnt = genPnt, metric = metric, dist = dist
  ))

}

##### other things #############################################################

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
