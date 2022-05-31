### math functions


##### model spaces: Euclidean, Sphere, and Poincare disk #######################

# what is needed:
#   coordinate system and random generator
#   corresponding metric tensor and geodesic distance
#   embedding into Euclidean space
#   (maybe) exponential/logarithmic map

#' Euclidean space
#'
#' @param d dimension
#'
#' @return a list
#' @export
#'
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

#' d-dimensional hyperbolic space
#'
#' @param d dimension
#' @param r radius
#' @param model \code{"ball"} or \code{"half-space"}
#'
#' @return a list
#' @details See Theorem 3.7, problem 6-3 of Lee (2018), and perhaps, the
#' wiki page of Poincare half-plane model.
#'
#' @export
#'
#' @examples
#' d <- 2
#' manifold <- spaceHyperbolic(d = d, 'ball')
#' # manifold <- spaceHyperbolic(d = d, 'half')
#' df.ell <- with(manifold, {
#'   getMetricEllipsoid(genGrid(11), metric, radius = 0.075)
#' })
#' head(df.ell)
#' ggplot2::ggplot(df.ell, ggplot2::aes(x, y, group = idx.pnt)) +
#'   ggplot2::geom_path() + ggplot2::coord_fixed() +
#'   ggplot2::geom_path(
#'     data = as.data.frame(
#'       getGeodesic(manifold$metric, d = d)(c(0.5, 0.5), c(-0.5, 0.05))
#'     ),
#'     mapping = ggplot2::aes(x1, x2), color = 'red', inherit.aes = F
#'   )
spaceHyperbolic <- function(d, r = 1, model = 'ball'){

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

#' d-dimensional sphere
#'
#' @param d dimension
#' @param r radius
#' @param chart \code{"stereographic"} or \code{"polar"}
#'   the coordinate chart to use
#'
#' @return a list
#' @details For stereographic projection, north pole \eqn{(0, ..., 0, 1} is
#' excluded, c.f. problem 1-7, Lee (2013), and pp. 58 - 61, Lee (2018).
#' \cr
#' For polar coordinate, the embedding of
#' \eqn{(r, \theta_1, \dots, \theta_d) \in S^d} to
#' \eqn{(x_1, \dots, x_{d+1}) \in R^{d+1}} follows
#' \eqn{x_1 = r \cos \theta_1}, \eqn{x_{d + 1} = r \prod_{j=1}^d \sin \theta_j},
#' and \deqn{
#'   x_i = r \cos \theta_i \prod_{j=1}^{i-1}\sin \theta_j, i = 2, \dots, d - 1.
#' }
#'
#' @examples
#' d <- 2
#' manifold <- spaceSphere(d = d, chart = 'polar')
#' # manifold <- spaceSphere(d = d, chart = 'stere')
#' df.ell <- with(manifold, {
#'   getMetricEllipsoid(genGrid(11), metric, radius = 0.075)
#' })
#' head(df.ell)
#' ggplot2::ggplot(df.ell, ggplot2::aes(x, y, group = idx.pnt)) +
#'   ggplot2::geom_path() + ggplot2::coord_fixed() +
#'   ggplot2::geom_path(
#'     data = as.data.frame(
#'       getGeodesic(manifold$metric, d = d)(c(2.5, 0), c(0, 5))
#'     ),
#'     mapping = ggplot2::aes(x1, x2), color = 'red', inherit.aes = F
#'   ) +
#'   ggplot2::labs(subtitle = 'x = 0 and x = pi are north and south pole.')
spaceSphere <- function(d, r = 1, chart = 'stereographic'){

  # d-dimensional sphere under ??? coordinates
  # coord chart: polar or stereographic projection from North pole
  # c.f. pp. 58 - 61 of Lee (2018) for metric (especially eq (3.7))
  # (also https://math.stackexchange.com/questions/402102/what-is-the-metric-tensor-on-the-n-sphere-hypersphere)

  stopifnot(d >= 1)
  chart <- match.arg(chart, c('stereographic', 'polar'))

  if(chart == 'stereographic'){

    genPnt <- function(n){
      # random point generator, returns a matrix, one row = one point
      stopifnot(n > 0)
      return(3 * r * matrix(runif(n * d, -1, 1), ncol = d, nrow = n))
    }

    genGrid <- function(n){
      # give a grid with equally spaced (in coordinates) n^d points
      stopifnot(n > 0)
      grid.margin <- 2 * r * seq(-1, 1, length.out = n)
      grid <- do.call(expand.grid, rep(list(grid.margin), d))
      grid <- as.matrix(grid)
      colnames(grid) <- NULL
      return(grid)
    }

    # metric tensor
    metric <- function(x){
      # x in stereographic projection coord
      return(
        4 * r ^ 4 * diag(d) / (sum(x ^ 2) + r^2) ^ 2
      )

    }

    # geodesic distance
    dist <- function(x, y){
      # function to compute geodesic distance on sphere space
      # x, y in stereographic projection coord, one row for one point

      if(!is.matrix(x)) x <- matrix(x, nrow = 1)
      if(!is.matrix(y)) y <- matrix(y, nrow = 1)
      stopifnot(all(dim(x) == dim(y)))
      stopifnot(ncol(x) == d)

      # embed to Euclidean space R^{d+1} and compute
      cart.x <- stereo2cart(x)
      cart.y <- stereo2cart(y)
      return(
        r * acos(rowSums(cart.x * cart.y) / r ^ 2)
      )
    }

  }else{

    # polar coordinates

    genPnt <- function(n){
      # random point generator, returns a matrix, one row = one point
      # here we let the thetas to be inside (0, pi)
      stopifnot(n > 0)
      return(pi * matrix(runif(n * d), ncol = d, nrow = n))
    }

    genGrid <- function(n){
      # give a grid with equally spaced (in coordinates) n^d points
      stopifnot(n > 0)
      grid.margin <- pi * seq(.Machine$double.eps ^ (1/2), 1, length.out = n)
      grid <- do.call(expand.grid, rep(list(grid.margin), d))
      grid <- as.matrix(grid)
      colnames(grid) <- NULL
      return(grid)
    }

    # metric tensor
    metric <- function(x){
      # x in polar coord
      return(
        r ^ 2 * diag( c(1, cumprod(sin(x) ^ 2)[-d] ) )
      )
    }

    # geodesic distance
    dist <- function(x, y){
      # function to compute geodesic distance on sphere space
      # x, y in polar coord, one row for one point

      if(!is.matrix(x)) x <- matrix(x, nrow = 1)
      if(!is.matrix(y)) y <- matrix(y, nrow = 1)
      stopifnot(all(dim(x) == dim(y)))
      stopifnot(ncol(x) == d)

      # embed to Euclidean space R^{d+1} and compute
      cart.x <- pol2cart(theta = x, r = r)
      cart.y <- pol2cart(theta = y, r = r)
      return(
        r * acos(rowSums(cart.x * cart.y) / r ^ 2)
      )
    }

  }


  return(list(
    genPnt = genPnt, genGrid = genGrid, metric = metric, dist = dist
  ))

}

pol2cart <- function(x, theta, r){

  # transform polar coord -> cartesian in R^n, n >= 2.
  # args:
  # x = a matrix whose row is (r, theta1, ...)
  # or theta and r separately
  # returns:
  # a matrix whose row is (x1, x2, ...)
  # x1 = r * cos(theta1),
  # xi = r * cos(theta_i) * Prod_{1 <= j <= i - 1} sin(theta_j),
  # xn = r * Prod_{1 <= j <= n - 1} sin(theta_j).
  # for i = 2, ..., n - 1.

  if(!missing(x)){
    if (!is.matrix(x)) x <- matrix(x, nrow = 1)
    r <- x[, 1]
    theta <- x[, -1, drop = FALSE]
  }else{
    r <- r
    if (!is.matrix(theta)) theta <- matrix(theta, nrow = 1)
    if(length(r) == 1)
      r <- rep(r, nrow(theta))
    else
      stopifnot(length(r) == nrow(theta))
  }

  n.dim <- ncol(theta) + 1
  stopifnot(n.dim >= 2)

  n.pts <- length(r)
  res <- matrix(0, ncol = n.dim, nrow = n.pts)
  tm <- r
  for (i in seq_len(n.dim - 1)) {
    res[, i] <- tm * cos(theta[, i])
    tm <- tm * sin(theta[, i])
  }
  res[, n.dim] <- tm
  return(res)
}

stereo2cart <- function(x){
  # translate stereographic projection coord to Cartisian coord
  # c.f. problem 1-7 of Lee (2013)
  # x: an array or matrix w/ 1row = 1point

  if(!is.matrix(x)) x <- matrix(x, nrow = 1)
  normsq.x <- rowSums(x ^ 2)

  return(
    cbind(2 * x, normsq.x - 1) / (normsq.x + 1)
  )
}
# stereo2cart(rep(0, 2))
# stereo2cart(rep(1, 2))
cart2stereo <- function(x){
  # computes stereographic projection from north pole (0, ..., 0, 1)
  if(!is.matrix(x)) x <- matrix(x, nrow = 1)
  d <- ncol(x)
  stopifnot(d >= 2)
  return(
    x[, -d, drop = F] / (1 - x[, d])
  )
}
# stereo2cart(rep(0, 2)) %>% cart2stereo
# stereo2cart(rep(1, 2)) %>% cart2stereo

# likely no need
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
getGeodesic <- function(metric, christ, d){
  # get a function computes geodesic

  stopifnot((!missing(metric) | !missing(christ)) & !missing(d))
  stopifnot(d >= 1)

  if(missing(christ)){

    stopifnot(inherits(metric, 'function'))
    geo.eq <- getGeoEq(metric = metric, d = d)

  }else{

    stopifnot(inherits(christ, 'function'))
    geo.eq <- getGeoEq(christ = christ, d = d)

  }


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
#' @param christ christoffel symbol (function)
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
getGeoEq <- function(metric, christ, d){

  stopifnot((!missing(metric) | !missing(christ)) & !missing(d))

  if(missing(christ)){
    metric <- metric
    d <- d

    stopifnot(inherits(metric, 'function'))
    stopifnot(d >= 1)

    christ.f <- getChristoffel(metric, d)
  }else{

    christ.f <- christ

  }

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
#' christ.f <- getChristoffel(manifold$metric, d = 2)
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
      aperm(arr.jacob, c(3, 1, 2)) + aperm(arr.jacob, c(1, 3, 2)) - arr.jacob

    # augment to arr.diff[i, j, k, l] = arr.diff[i, j, l] for all k = 1, ..., d
    aug.diff <- array(arr.diff, dim = rep(d, 4))
    aug.diff <- aperm(aug.diff, c(1, 2, 4, 3))
    # # check correctness
    # apply(aug.diff, 3, function(tm) identical(arr.diff, tm)) %>%
    #   all %>% stopifnot

    # augment inverse metric matrix to aug.inv.met[i, j, k, l]=inv.metric[k, l]
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

# d <- 2
# metric <- spaceSphere(d = d)$metric
# get04Curvature(metric, d)(rep(1, d))

get04Curvature <- function(metric, d) {

  # the (0,4) curvature tensor, c.f. eq. (7.7) & (7.8) of Lee (2018)

  metric <- metric
  d <- d

  c13.f <- get13Curvature(metric, d)

  res.f <- function(x) {
    # browser()
    # the (0,4) curvature tensor, c.f. eq. (7.7) & (7.8) of Lee (2018)

    # cvt.13[i,j,k,l] = R_{i,j,k}^l
    cvt.13 <- c13.f(x)
    # aug.cvt.13[i,j,k,l,m] = cvt.13[i,j,k,l]
    aug.cvt.13 <- array(cvt.13, rep(d, 5))
    # permute aug.cvt.13[i,j,k,l,m] <- aug.cvt.13[i,j,k,m,l] = cvt.13[i,j,k,m]
    # aug.cvt.13 <- aperm(aug.cvt.13, c(1, 2, 3, 5, 4))
    aug.cvt.13 <- reindex(aug.cvt.13, to = 'ijklm', from = 'ijkml')

    # aug.met[i,j,k,l,m] = metric[i,j]
    aug.met <- array(metric(x), rep(d, 5))
    # reshape2::melt(aug.met)
    # permute, aug.met[i,j,k,l,m] <- aug.met[l,m,i,j,k] = metric[l,m]
    # aug.met <- aperm(aug.met, perm = c(3, 4, 5, 1, 2))
    aug.met <- reindex(aug.met, to = 'ijklm', from = 'lmijk')
    # aug.met[,,,1,1] #?????
    # reshape2::melt(aug.met)
    # tst <- array(seq(32), rep(2, 5));
    # reshape2::melt(tst) %>% select(Var4, Var5, Var1, Var2, Var3, value) %>%
    #   arrange(Var3, Var2, Var1, Var5, Var4)
    # reshape2::melt(aperm(tst, perm = c(4,5,1,2,3)))

    # compute (0,4)-tensor
    res.arr <- rowSums(
      # afterProd[i,j,k,l,m] = aug.met[i,j,k,l,m]*aug.cvt.13[i,j,k,l,m] =
      #..=metric[l,m] * cvt.13[i,j,k,m]
      aug.met * aug.cvt.13,
      dims = 4
    )

    return(res.arr)

  }

  return(res.f)

}

get13Curvature <- function(metric, d){
  # compute the (1,3)-curvature tensor

  stopifnot(!missing(metric) & !missing(d))

  metric <- metric
  d <- d

  stopifnot(inherits(metric, 'function'))
  stopifnot(d >= 1)

  christ <- getChristoffel(metric, d)
  # nloptr::nl.jacobian(rep(1, d), christ) # [i, j] will be \partial_j f_i

  res.f <- function(x) {

    # compute the (1,3)-curvature tensor at x, following eq.(7.4) of Lee (2018)

    # x <- rep(1, d)

    # partial.christ[i, j] = \partial_j of vectorized_christ[i]
    partial.christ <- nloptr::nl.jacobian(x, christ)
    # into array, hence partial.christ[i, j, k, l] = \partial_l Christ[i, j, k]
    partial.christ <- array(partial.christ, dim = rep(d, 4))
    # rearrange into [i, j, k, l] = \partial_i Christoffel[j, k, l]
    # partial.christ <- aperm(partial.christ, c(4, 1, 2, 3))
    partial.christ <- reindex(partial.christ, to = 'ijkl', from = 'jkli')

    # now the difference of
    #..\partial_i Christ[j, k, l] - \partial_j Christ[i, k, l]
    arr.diff <- partial.christ -
      # after_permute[i,j,k,l]=partial.christ[j,i,k,l]
      # aperm(partial.christ, c(2, 1, 3, 4))
      reindex(partial.christ, to = 'ijkl', from = 'jikl')

    # augment christ [i, j, k, l, m] = Christ[i, j, k]
    aug.christ <- array(christ(x), rep(d, 5))
    # product term [i, j, k, l, m] = Christ[j, k, m] * Christ[i, m, l]
    arr.prod <-
      # after_permute[i,j,k,l,m]=aug.christ[j,k,m,i,l]=Christ[j,k,m]
      # aperm(aug.christ, c(2, 3, 5, 1, 4)) *
      reindex(aug.christ, to = 'ijklm', from = 'jkmil') *
      # after_permute[i,j,k,l,m]=aug.christ[i,m,l,j,k]=Christ[i,m,l]
      # aperm(aug.christ, c(1, 5, 4, 2, 3))
      reindex(aug.christ, to = 'ijklm', from = 'imljk')
    # sum over m, so
    #..arr.prod[i, j, k, l] = sum_m ( Christ[j, k, m] * Christ[i, m, l] )
    arr.prod <- rowSums(arr.prod, dims = 4)

    # the (1,3)-curvature tensor
    res.arr <- arr.diff + arr.prod -
      # after_permute[i,j,k,l]=arr.prod[j,i,k,l]
      reindex(arr.prod, to = 'ijkl', from = 'jikl')
      # aperm(arr.prod, c(2, 1, 3, 4))

    return(res.arr)

  }

  return(res.f)

}

getSecCurvature <- function(metric, riem.cvt, d){

  # computes sectional curvature, c.f. pp.250 - 255, Lee (2018)
  # metric & riem.cvt: metric anr (0, 4)-curvature tensor,
  # must be both function or both tensor array. riem.cvt can be missing if
  # metric is a function.
  # d: dimension
  # returns a function (target, v, w) computing sectional curvature
  # at target of the plane spanned by v and w.

  # if(inherits(riem.cvt, 'function')){
  #   riem.cvt <- riem.cvt
  #   res.f <- function(target, ...) getSecCurvature(riem.cvt(target), ...)
  #   return(res.f)
  # }

  if(missing(riem.cvt)) riem.cvt <- NULL

  if(!inherits(metric, 'function') & !inherits(riem.cvt, 'function')){

    mat.met <- metric
    tsr <- riem.cvt

    # if input tensor array
    res.f <- function(v, w){

      stopifnot(length(v) == d)
      stopifnot(length(v) == length(w))

      stopifnot(all(dim(mat.met) == rep(length(v), 2)))
      stopifnot(all(
        dim(tsr) == rep(length(v), 4)
      ))

      # Riemann (0, 4)-curvature R(v, w, w, v)
      tm.nu <- getTensorAction(tsr, v, w, w, v)
      # wedge norm squared
      tm.de <-
        getTensorAction(mat.met, v, v) * getTensorAction(mat.met, w, w) -
        getTensorAction(mat.met, v, w) ^ 2

      return(tm.nu / tm.de)

    }

  }else{

    stopifnot(inherits(metric, 'function'))
    # riem.cvt <- get04Curvature(metric, d)
    if(is.null(riem.cvt)) riem.cvt <- get04Curvature(metric, d)

    res.f <- function(target, v, w) getSecCurvature(
      metric = metric(target), riem.cvt = riem.cvt(target), d = d
    )(v, w)

  }
  return(res.f)


  # res.f <- function(target, v, w){
  #
  #   stopifnot(length(v) == d)
  #   stopifnot(length(v) == length(w))
  #
  #   mat.met <- metric(target)
  #   tsr <- riem.cvt(target)
  #   stopifnot(all(dim(mat.met) == rep(length(v), 2)))
  #   stopifnot(all(
  #     dim(tsr) == rep(length(v), 4)
  #   ))
  #
  #   # Riemann (0, 4)-curvature R(v, w, w, v)
  #   tm.nu <- getTensorAction(tsr, v, w, w, v)
  #   # wedge norm squared
  #   # idx <- vecQuadIdx(d = nrow(mat.met), idx.vec = T)
  #   # vec.v <- vecQuad(v)
  #   # vec.w <- vecQuad(w)
  #   tm.de <-
  #     getTensorAction(mat.met, v, v) * getTensorAction(mat.met, w, w) -
  #     getTensorAction(mat.met, v, w) ^ 2
  #
  #   return(tm.nu / tm.de)
  #
  # }


}

##### additional function ######################################################

# sigmoid function and its inverse
sigmoid.f <- function(x){
  1 / (1 + exp(-x))
}
sigmoid.f.inv <- function(x){
  log(x / (1 - x))
}
