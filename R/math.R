### math functions


##### model spaces: Euclidean, Sphere, and Poincare disk #######################

# what is needed:
#   coordinate system and random generator
#   corresponding metric tensor and geodesic distance
#   embedding into Euclidean space
#   (maybe) exponential/logarithmic map
#' Title
#'
#' @param d
#'
#' @return
#' @export
#'
#' @examples
spaceEuclidean <- function(d){
  # d-dimensional Euclidean space under Cartesian coordinates

  stopifnot(d >= 1)

  genPnt <- function(n){
    # random point generator, returns a matrix, one row = one point
    stopifnot(n > 0)
    return(matrix(runif(n * d), ncol = d, nrow = n))
  }

  genGrid <- function(n){
    # give a grid with equally spaced (in coordinates) n^d points
    stopifnot(n > 0)
    grid.margin <- seq(0, 1, length.out = n)
    grid <- do.call(expand.grid, rep(list(grid.margin), d))
    grid <- as.matrix(grid)
    colnames(grid) <- NULL
    return(grid)
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
    genPnt = genPnt, genGrid = genGrid, metric = metric, dist = dist
  ))

}

#' Title
#'
#' @param d
#' @param model
#' @param r
#'
#' @return
#' @export
#'
#' @examples
spaceHyperbolic <- function(d, model, r = 1){

  ### d-dim Hyperbolic space with radius r, c.f. Theorem 3.7 of Lee (2018).

  model <- match.arg(model, c('ball', 'half-space'))

  if(model == 'ball'){

    ## Poincare ball model
    genPnt <- function(n){

      # random point generator, returns a matrix, one row = one point
      stopifnot(n > 0)
      vol.d.dim.ball <- pi^(d/2) / gamma(d / 2 + 1) * r^d
      n.fill <- ceiling(2 * n / vol.d.dim.ball * (2 * r)^d)
      tm <- matrix(runif(n.fill * d, -r, r), ncol = d, nrow = n.fill)
      tm <- tm[rowSums(tm^2) < r^2 - .Machine$double.eps ^ (1/2), , drop = F]
      tm <- head(tm, n)
      return(tm)

    }

    genGrid <- function(n){
      # give a grid with equally spaced (in coordinates) n^d points
      stopifnot(n > 0)
      grid.margin <- seq(-r, r, length.out = n)
      grid <- do.call(expand.grid, rep(list(grid.margin), d))
      tm <- as.matrix(grid)
      tm <- tm[rowSums(tm^2) < r^2 - .Machine$double.eps ^ (1/2), , drop = F]
      colnames(tm) <- NULL
      return(tm)
    }

    metric <- function(x){
      # x in Cartesian coord
      return(
        diag(
          rep(4 * r^4 / (r^2 - sum(x^2))^2, length(x))
        )
      )
    }
    dist <- function(x, y){

      # function to compute geodesic distance in Poincare ball model
      # x, y in Cartesian coord, one row for one point,
      # c.f. problem 6-3 of Lee (2018)

      if(!is.matrix(x)) x <- matrix(x, nrow = 1)
      if(!is.matrix(y)) y <- matrix(y, nrow = 1)
      stopifnot(all(dim(x) == dim(y)))
      rsq.x <- rowSums(x^2)
      rsq.y <- rowSums(y^2)
      stopifnot(all(rsq.x < r^2) & all(rsq.y < r^2))
      distsq.eucl <- rowSums((x - y)^2)

      return(
        r * acosh(
          1 + (2 * r^2 * distsq.eucl) / ((r^2 - rsq.x) * (r^2 - rsq.y))
        )
      )
    }

  }else{

    ## Poincare half-space model

    genPnt <- function(n){
      # random point generator, returns a matrix, one row = one point
      stopifnot(n > 0)
      return(matrix(runif(n * d), ncol = d, nrow = n))
    }

    genGrid <- function(n){
      # give a grid with equally spaced (in coordinates) n^d points
      stopifnot(n > 0)
      grid.margin <- seq(.Machine$double.eps ^ (1/3), 1, length.out = n)
      grid <- do.call(expand.grid, rep(list(grid.margin), d))
      grid <- as.matrix(grid)
      colnames(grid) <- NULL
      return(grid)
    }

    metric <- function(x){
      # x in Cartesian coord
      return(
        diag(
          rep(r^2 / x[length(x)]^2, length(x))
        )
      )
    }
    dist <- function(x, y){

      # function to compute geodesic distance in Poincare ball model
      # x, y in Cartesian coord, one row for one point,
      # c.f. problem 6-3 of Lee (2018) and Wiki: Poincare half-plane model

      if(!is.matrix(x)) x <- matrix(x, nrow = 1)
      if(!is.matrix(y)) y <- matrix(y, nrow = 1)
      stopifnot(all(dim(x) == dim(y)))

      d <- ncol(x)
      distsq.eucl <- rowSums((x - y)^2)
      prod.last.coord <- x[, d] * y[, d]

      return(
        r * acosh(
          1 + (distsq.eucl) / (2 * prod.last.coord)
        )
      )
    }

  }


  return(list(
    genPnt = genPnt, genGrid = genGrid, metric = metric, dist = dist
  ))


}

#' Title
#'
#' @param d
#'
#' @return
#' @export
#'
#' @examples
spaceHemisphere <- function(d){
  # d-dimensional hemisphere under ??? coordinates

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


locDist <- function(coord, dist, local.reach){
  # get some idea about magnitude of local distance

  coord.1 <- coord - local.reach / 2
  coord.2 <- coord + local.reach / 2

  return(dist(
    coord.1, coord.2
  ))

}

##### Riemannian geometry related functions ####################################

#' Compute geodesic curve
#'
#' @param metric metric function
#' @param d dimension
#'
#' @return a function
#'
#' @export
#'
#' @details The resulting function takes start point \code{start.pnt},
#' initial velocity \code{init.v}, an array of time \code{t},
#' then computes geodesic curve and return in a \code{1 + 2 * d}-col matrix,
#' where the columns are time, the curve, and the velocity.
#'
#' @examples
#' manifold <- spaceHyperbolic(d = 2, model = 'half')
#' geo.f <- getGeodesic(manifold$metric, d = 2)
#' # solve for geodesic curve from c(0.5, 0.5) w/ initial velocity c(1, 0)
#' geo <- geo.f(c(0.5, 0.5), c(1, 0))
#' plot(x = geo[, 'x1'], y = geo[, 'x2'], type = 'l')
getGeodesic <- function(metric, d){
  # get a function computes geodesic

  stopifnot(!missing(metric) & !missing(d))

  metric <- metric
  d <- d

  stopifnot(inherits(metric, 'function'))
  stopifnot(d >= 1)

  geo.eq <- getGeoEq(metric, d)
  res.f <- function(start.pnt, init.v, t = seq(0, 1, by = 0.01)){

    # geodesic curve from start w/ initial velocity, given grid of time
    stopifnot(length(start.pnt) == length(init.v))
    stopifnot(length(start.pnt) == d)
    stopifnot(all(diff(t) > 0))
    stopifnot(all(t >= 0))

    ode.res <- deSolve::ode(
      c(start.pnt, init.v),
      times = t, func = geo.eq
    )
    colnames(ode.res) <-
      c('time', sprintf('x%s', seq(d)), sprintf('v%s', seq(d)))
    return(ode.res)
  }

  return(res.f)

}

#' Construct ODE for geodesic equation
#'
#' @param metric metric function
#' @param d dimension
#'
#' @return a function representing geodesic equation following input metric that
#' can be passed to \code{deSolve::ode}.
#' @export
#'
#' @examples
#' manifold <- spaceHyperbolic(d = 2, model = 'half')
#' geo.eq <- getGeoEq(manifold$metric, d = 2)
#' geo.eq(0, c(rep(0.5, 2), c(1, 0)))
#' # solve for geodesic curve from c(0.5, 0.5) w/ initial velocity c(1, 0)
#' geo <- deSolve::ode(
#'   c(rep(0.5, 2), c(1, 0)), times = seq(0, 5, by = 0.01), func = geo.eq
#' )
#' plot(x = geo[, 2], y = geo[, 3], type = 'l')
getGeoEq <- function(metric, d){

  stopifnot(!missing(metric) & !missing(d))

  metric <- metric
  d <- d

  stopifnot(inherits(metric, 'function'))
  stopifnot(d >= 1)

  christ.f <- getChristoffel(metric, d)

  geo.eq <- function(t, y, ...) {

    # ODE for geodesic equation, t: time, y: vector of (x, v),
    # where x is current location, v is current velocity.

    vec.x <- y[seq(d)]
    vec.v <- y[d + seq(d)]

    arr.christ <- christ.f(vec.x)

    dx <- vec.v
    tm.arr <- array(outer(vec.v, vec.v), dim = dim(arr.christ))
    dv <- -1 * colSums( tm.arr * arr.christ, dims = 2 )

    return(list(c(dx, dv)))

  }

  return(geo.eq)

}


#' Compute Christoffel symbol given metric
#'
#' @param metric a metric function
#' @param d dimension
#'
#' @return a function with argument \code{x}.
#'
#' @export
#'
#' @details The resulting function computes Christoffel symbol and returns a
#' \code{d X d X d} array, whose \code{[i, j, k]} element is
#' \deqn{
#' \Gamma^k_{ij} = \frac{1}{2} g^{kl} (
#'   \partial_i g_{jl} + \partial_j g_{il} - \partial_l g_{ij}),
#' }
#' where \eqn{g^{kl}} is the \eqn{(k, l)}-element of the inverse of metric
#' tensor. See also Lee (2018) pp.123 - 124.
#' @examples
#' manifold <- spaceHyperbolic(d = 2, model = 'half')
#' christ.f <- get_Christoffel(manifold$metric, d = 2)
#' christ.f(rep(0.1, 2))
getChristoffel <- function(metric, d) {
  # numerically compute Christoffel symbol of the input metric

  stopifnot(!missing(metric) & !missing(d))

  metric <- metric
  d <- d

  stopifnot(inherits(metric, 'function'))
  stopifnot(d >= 1)

  christ.f <- function(x) {

    # computing Christoffel symbol, c.f. Lee(2018) pp.124-124.
    # returns an array whose [i, j, k] element is \Gamma^k_{ij}

    # x <- rep(0.1, d)

    # inverse metric tensor
    inv.metric <- MASS::ginv(metric(x))

    mat.jacob <- nloptr::nl.jacobian(x, metric)
    # now mat.jacob[i, j] = \partial_j g_{pq} w/ i = p = (q-1)*d
    # re-arrange into d X d X d array w/ arr.jacob[i, j, l] = \partial_l g_{ij}
    arr.jacob <- array(mat.jacob, dim = rep(d, 3))

    # compute arr.diff[i, j, l] =
    #     \partial_i g_{jl} + \partial_j g_{il} - \partial_l g_{ij}
    arr.diff <-
      aperm(arr.jacob, c(2, 3, 1)) + aperm(arr.jacob, c(3, 2, 1)) - arr.jacob

    # augment to arr.diff[i, j, k, l] = arr.diff[i, j, l] for all k = 1, ..., d
    aug.diff <- array(arr.diff, dim = rep(d, 4))
    aug.diff <- aperm(aug.diff, c(1, 2, 4, 3))
    # # check correctness
    # apply(aug.diff, 3, function(tm) identical(arr.diff, tm)) %>%
    #   all %>% stopifnot

    # augment inverse metric matrix to aug.inv.met[i, j, k, l] = inv.metric[k, l]
    aug.inv.met <- array(inv.metric, dim = rep(d, 4))
    aug.inv.met <- aperm(aug.inv.met, c(3, 4, 1, 2))
    # # check correctness
    # apply(aug.inv.met, c(1, 2), function(tm) {
    #   identical(inv.metric, tm)
    # }) %>% all %>% stopifnot

    # the Christoffel symbol arr.christ[i, j, k] = \Gamma^k_{ij}
    arr.christ <- rowSums(aug.inv.met * aug.diff / 2, dims = 3)
    # arr.christ[,, 2]
    # arr.christ[,, 1]

    return(arr.christ)

  }

  return(christ.f)

}

# No need now
# derivation <- function(f){
#   # compute derivation functions (in coord representation) of f
#   # basically gives you \partial_1 f, \dots, \partial_d f, in function form
#
#   # wrapper over nloptr::nl.grad / nloptr::nl.jacobian
#
#   res.f <- function(x) nloptr::nl.grad(x, f)
#   return(res.f)
#
# }

##### vectorization of quadratic forms #########################################

vecQuadIdx <- function(d, idx.vec = F){
  # generate combinations of seq(d) taken 2 at a time
  # with replacement.
  # Returns: a 2 X d(d-1)/2 where each column is a
  # combination.
  stopifnot(d >= 2)
  res <- cbind(
    # diag elements
    t(matrix(seq(d), nrow = d, ncol = 2)),
    utils::combn(d, 2) # non-diag elements
  )

  if(!idx.vec){
    return(res)
  }else{
    return(
      res[2, ] + d * (res[1, ] - 1)
    )
  }

}

symMat2Vec <- function(mat) {#vecSymMat(mat)
# vecSymMat <- function(mat){
  # vectorize a symmetric matrix following vecQuadIdx
  # input a d-by-d symmetric matrix
  # output a d * (d + 1) / 2 array

  stopifnot(isSymmetric(mat))
  d <- nrow(mat)
  idx <- vecQuadIdx(d)
  idx <- idx[2, ] + d * (idx[1, ] - 1)
  return(as.numeric(mat[idx]))
}
# symMat2Vec(diag(3))
# symMat2Vec(matrix(seq(9), 3) + t(matrix(seq(9), 3)))

vec2SymMat <- function(vec){
  # construct symmetric matrix following vecQuadIdx
  d <- (sqrt(8 * length(vec) + 1) - 1) / 2
  mat <- matrix(NA, d, d)
  mat[vecQuadIdx(d, T)] <- vec
  mat[is.na(mat)] <- mat[is.na(t(mat))]
  # mat <- (mat + t(mat)) / 2
  # diag(mat) <- vec[seq(d)]
  return(mat)
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
#' @details compute \code{t(mat.x) \%*\% mat.a \%*\% mat.x}, so the (i, j)
#' element is \code{mat.x[, i] \%*\% mat.a \%*\% mat.x[, j]}.
#' @export
#' @md
#' @examples TBD
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
