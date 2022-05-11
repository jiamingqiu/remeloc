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
      return(3 * r * matrix(runif(n * d), ncol = d, nrow = n))
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

##### vectorization of Riemann curvature tensor ################################

# getTensorAction <- function(tensor){
#   # vectorize x for later action from Riemann curvature (0,4)-tensor
#   # TBD
#   res.f <- function(...){
#     as.name(match.call())
#   }
# }

modelMatRiemCvt <- function(df.idx, x, y, z, w){

  # construct design matrix for local regression on Riemann (0,4)-curvature
  # df.idx: result from vecRiemCvtIdx
  # ...: 4 tensor arrays, no sanity check, make sure dimension match!

  # ls.args <- as.list(match.call())
  # ls.args <- ls.args[-1]
  # ls.args <- ls.args[!(names(ls.args) %in% 'df.idx')]
  # stopifnot(length(ls.args) == 4)

  arr.r <- with(df.idx, {
    x[, i1] * y[, i2] * z[, i3] * w[, i4]
    # ls.args[[1]][, i1] * ls.args[[2]][, i2] *
    #   ls.args[[3]][, i3] * ls.args[[1]][[4]][, i4]
  }) # arr.r[, j] = term in df.idx[j, ]
  # apply(arr.r, 2, sd) %>% log10
  vec.cvt.basis <- as.matrix(
    df.idx[, stringr::str_detect(names(df.idx), 'e\\d+'), drop = F]
  )
  # put the coef for basis into columns
  mat.r <- -1 / 3 * apply(vec.cvt.basis, 2, function(x) {
    rowSums(arr.r * x[col(arr.r)])
  })

  return(mat.r)

}

#' Index table for vectorizing algebraic curvature tensor
#'
#' @param d dimension of vector space
#'
#' @return a data frame
#' @export
#'
#' @details See \link{\code{vecTensorIdx}} for details. See p212, Lee (2018)
#' for definition.
#'
#' @examples
#' d <- 3
#' manifold <- spaceHyperbolic(d, model = 'half')
#' cvt04.f <- get04Curvature(manifold$metric, d = d)
#' pnt <- manifold$genPnt(10)
#' true.cvt <- apply(pnt, 1, cvt04.f) # flatten tensors
#' df.idx <- vecRiemCvtIdx(d, drop.zero = F)
#' # coefficients, the basis are orthonormal
#' mat.basis <- as.matrix(
#'   df.idx[, stringr::str_detect(names(df.idx), 'e\\d+'), drop = F]
#' )
#' coef.cvt <- apply(true.cvt, 2, function(x){
#'   colSums(x * mat.basis)
#' })
#' coef.cvt <- matrix(coef.cvt, nrow = ncol(mat.basis))
#' reconstruct.cvt <- mat.basis %*% coef.cvt
#' all.equal(reconstruct.cvt, true.cvt)
#' see.cvt <- cvt04.f(pnt[1, ])
#' for(idx in seq(nrow(df.idx))){
#'   stopifnot(all.equal(
#'     reconstruct.cvt[idx, 1],
#'     with(
#'       df.idx[idx, sprintf('i%s', seq(4))],
#'       see.cvt[i1, i2, i3, i4]
#'     )
#'   ))
#' }
vecRiemCvtIdx <- function(d, ...) {
  # generate index of d-dim Riemann (0,4)-Curvature tensor following symmetries

  sym.eq <- c(
    'ijkl = -jikl' # means T(w, x, y, z) = T(x, w, y, z)
    , 'ijkl = - ijlk'
    , 'ijkl = klij'
    , 'ijkl = -jkil - kijl'
  )
  res.dim <- d^2 * ( d^2 - 1 ) / 12
  df.subscript <- vecTensorIdx(d = d, k = 4, sym.eq = sym.eq, ...)
  stopifnot(ncol(df.subscript) == res.dim + 4 + 1 + 1)

  return(df.subscript)

}

#' Index table for vectorizing tensor/multidimensional array w.r.t. symmetries
#'
#' @param d dimension of vector space
#' @param k dimension of indices
#' @param sym.eq equalities defining symmetries
#' @param drop.zero whether to drop zeros in return
#'
#' @return a data frame
#' @export
#'
#' @details The each row in the returned data frame represents one element in
#' the desired tensor, and with the following columns:
#' \describe{
#'   \item{\code{idx}}{The index of the element if flatten.}
#'   \item{\code{i1, ..., ik}, and \code{subscript}}{Subscripts.}
#'   \item{\code{e1, ...}}{Basis used for representing.}
#' }
#' See \code{\link{decodeSymEq}} for how to specify symmetries.
#'
#' @examples
#' # the algebraic curvature tensor
#' vecTensorIdx(2, 4, sym.eq = c(
#'   'ijkl = -jikl' , 'ijkl = - ijlk',
#'   'ijkl = klij', 'ijkl = -jkil - kijl'
#' ))
vecTensorIdx <- function(d, k, sym.eq, drop.zero = T) {

  # generate index for vectorization of tensor

  # d <- 2
  # k <- 4
  # sym.eq <- list(
  #   'ijkl = -jikl' # means T(w, x, y, z) = T(x, w, y, z)
  #   , 'ijkl = - ijlk'
  #   , 'ijkl = klij'
  #   , 'ijkl = -jkil - kijl'
  # )
  # # res.dim <- d^2 * ( d^2 - 1 ) / 12

  stopifnot(d >= 1)
  stopifnot(k >= 1)

  df.subscript <- expand.grid(rep(list(seq(d)), k))
  names(df.subscript) <- sprintf('i%s', seq(k))
  mat.subscript <- as.matrix(df.subscript)
  df.subscript$subscript <- apply(mat.subscript, 1, function(idx){
    paste(idx, collapse = ',')
  })

  if(missing(sym.eq)){
    # if no symmetries, use all elements
    vec.basis <- diag(nrow(df.subscript))
  }else{
    # if(!is.list(sym.eq)) sym.eq <- list(sym.eq)

    # list of decoder according to rules
    ls.rule <- lapply(sym.eq, decodeSymEq)
    # constructing linear system per subscript combination accordingly
    ls.lin.sys <- lapply(ls.rule, function(rule) {
      apply(mat.subscript, 1, function(x) {
        rule(as.numeric(x), lin.sys = T)
      })
    })
    ls.lin.sys <- unlist(ls.lin.sys, recursive = F)
    # translate to coefficient vector
    ls.coef <- lapply(ls.lin.sys, function(lin.sys){
      where <- match(
        sapply(lin.sys$idx, paste, collapse = ','),
        df.subscript$subscript
      )
      arr.coef <- rep(0, nrow(df.subscript))
      arr.coef[where] <- lin.sys$coef
      arr.coef
    })
    # mat.coef is the coef matrix corr. sym.eqs, <==> mat.coef %*% x = 0
    mat.coef <- do.call(rbind, ls.coef)
    # use basis of its null space
    vec.basis <- pracma::null(mat.coef)
  }
  # some trimming due to rounding error
  vec.basis[abs(vec.basis) < .Machine$double.eps^(1/2)] <- 0
  # maybe no need
  # vec.basis <- apply(vec.basis, 2, function(x) x / max(abs(x)))
  # vec.basis %>% print(digits = 2)

  # combine and return
  df.basis <- as.data.frame(vec.basis)
  names(df.basis) <- sprintf('e%s', seq(ncol(df.basis)))
  df.subscript <- cbind(
    idx = seq(nrow(df.subscript)),
    df.subscript,
    df.basis
  )
  if(drop.zero){
    # only keep non-zero
    idx.keep <- rowSums(vec.basis != 0) >= 1

    return(
      df.subscript[idx.keep, , drop = F]
    )
  }else{
    return(df.subscript)
  }

}

#' Decoding symmetry equalities for tensor.
#'
#' @param eq equality determining symmetries
#'
#' @return a function
#'
#' @details
#' The equality specifies symmetries of tensor in terms of subscripts.
#' For example, given a (0,4)-tensor \eqn{R}, setting \code{eq} being
#' \code{'ijkl = -jikl'} means
#' \eqn{R_{ijkl} = -R_{jikl}}.
#' Similarly, \code{'ijkl = -jkil - kijl'} means
#' \eqn{R_{ijkl} = -R_{jkil} - R_{kijl}} (algebraic Bianchi identity).
#' There should be only one term to left hand side of \code{eq}, while the terms
#' to the right hand side should use identical indices (i,j,k,l, per say).
#' Also, indices for different dimension should be unique.
#' Currently only addition/subtraction is recognized.
#' \cr
#' The returned function takes argument \code{idx}, \code{name},
#' and \code{lin.sys}. See following examples for details.
#'
#' @export
#'
#' @examples
#' x <- array(seq(256), rep(4, 4))
#' # algebraic Bianchi identity for (0,4)-curvature tensor
#' res.f <- decodeSymEq('ijkl = -jkil - kijl')
#' x[1, 2, 3, 4]
#' res.f(c(1, 2, 3, 4))
#' x[2, 1, 3, 4]
#' - x[2, 3, 1, 4] - x[3, 1, 2, 4]
#' x[1, 2, 3, 4]
#' eval(parse(text = res.f(c(1, 2, 3, 4))))
#' x[1, 2, 3, 4]
#' res.f(c(1, 2, 3, 4), lin.sys = T)
decodeSymEq <- function(eq) {

  # decode symmetric using input subscript equality
  # gives a function that returns string "x[1, 2, 3, 4] <- x[2, 1, 3, 4]"

  # eq <- 'ijkl = - jikl'
  # eq <- 'ijkl = -jkil - kijl'
  eq <- stringr::str_remove_all(eq, ' ')
  eq <- stringr::str_split(eq, '=')[[1]]

  eq.lhs <- eq[1]
  stopifnot(!stringr::str_detect(eq.lhs, '(\\+)|(\\-)'))
  eq.rhs <- eq[2]

  # stringr::str_split('jklm+abcd-klmv', '.{0}(?=(\\+)|(\\-))') # magic...
  eq.rhs.terms <- stringr::str_split(eq.rhs, '.{0}(?=(\\+)|(\\-))')[[1]]
  eq.rhs.terms <- eq.rhs.terms[sapply(eq.rhs.terms, stringr::str_length) > 0]
  eq.rhs.sign <- stringr::str_extract(eq.rhs.terms, '(\\+)|(\\-)')
  eq.rhs.sign[is.na(eq.rhs.sign)] <- ''
  eq.rhs.terms <- stringr::str_remove_all(eq.rhs.terms, '(\\+)|(\\-)')

  eq.rhs.scr <- sapply(eq.rhs.terms, function(x) {
    stringr::str_split(x, '')[[1]]
  }, simplify = F)
  eq.lhs.scr <- stringr::str_split(eq.lhs, '')[[1]]

  # check index consistency
  stopifnot(!any(duplicated(eq.lhs.scr))) # unique subscript for each dim
  stopifnot(all(
    sapply(eq.rhs.scr, function(x) identical(
      stringr::str_sort(x), stringr::str_sort(eq.lhs.scr)
    ))
  ))

  eq.rhs.idx <- lapply(eq.rhs.scr, function(x) match(x, eq.lhs.scr))
  # construct resulting "deparse" function for lhs
  # eq.lhs
  # eq.rhs.sign
  # eq.rhs.idx
  res.f <- function(idx, name = 'x', lin.sys = F) {

    # idx: integer vector for index, 1, 2, 3, ...
    # name: the name of the array to handle
    # lin.sys: if full.assign = F, to return coef used for linear system or not
    # legacy: full.assign: to return "x[1,2] <- x[2,1]" or just the RHS "x[2,1]"


    stopifnot(length(idx) == length(eq.lhs.scr))
    stopifnot(is.character(name))

    # repermute according to prespecified rules
    ls.rhs.idx <- lapply(eq.rhs.idx, function(x) idx[x])
    rhs.str <- lapply(ls.rhs.idx, function(x){
      sprintf('%s[%s]', name, paste(x, collapse = ', '))
    })
    # add sign
    rhs.str <- mapply(
      function(term, sign) sprintf('%s %s', sign, term),
      term = rhs.str, sign = eq.rhs.sign
    )
    # combine
    rhs.str <- paste(rhs.str, collapse = '')
    if(
      length(eq.rhs.sign) == 1 & eq.rhs.sign[1] == '-' &
      all(ls.rhs.idx[[1]] == idx)
    ) {
      # if x[1,2,3,3] <- - x[1,2,3,3]
      rhs.str <- '0'
    }

    # if(!full.assign){
      if(lin.sys){
        ls.idx <- c(list(idx), ls.rhs.idx)
        names(ls.idx) <- c(eq.lhs, names(ls.rhs.idx))
        coef <- ifelse(eq.rhs.sign == '-', -1, 1)
        coef <- c(1, -1 * coef)
        # browser();QWER
        # use unique index only
        unique.idx <- unique(ls.idx)
        tbl <- match(as.character(ls.idx), as.character(unique.idx))
        # update
        coef <- as.numeric(sapply(split(coef, tbl), sum))
        return(list(coef = coef, idx = unique.idx))
      }#else
    #     return(rhs.str)
    # }

    lhs.str <- sprintf('%s[%s]', name, paste(idx, collapse = ', '))
    return(sprintf("%s <- %s", lhs.str, rhs.str))

  }

  return(res.f)

}

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

##### additional function ######################################################

# sigmoid function and its inverse
sigmoid.f <- function(x){
  1 / (1 + exp(-x))
}
sigmoid.f.inv <- function(x){
  log(x / (1 - x))
}

prodKN <- function(h, k){

  # Kulkarni--Nomizu product, c.f. eq (7.37) of Lee (2018)
  # not tested, use with care

  if(missing(k)) k <- h

  stopifnot(identical(class(h), class(k)))
  if(inherits(h, 'function')){
    res.f <- function(...){
      prodKN(h(...), k(...))
    }
    return(res.f)
  }

  stopifnot(all(dim(h) == dim(k)))
  d <- dim(h)[1]
  stopifnot(all(dim(h) == d))
  stopifnot(length(dim(h)) == 2)

  # augment h and k
  aug.h <- array(h, dim = rep(d, 4)) # aug.h[i,j,l,m]=h[i,j]
  aug.k <- array(k, dim = rep(d, 4)) # aug.k[i,j,l,m]=k[i,j]

  # arr.prod[i,j,l,m] = h[i,m] * k[j,l]
  arr.prod <-
    #..[i,j,l,m] = aug.h[i,m,j,l], <--> aug.h[i,j,l,m] = ..[i,l,m,j]
    aperm(aug.h, c(1, 3, 4, 2)) *
    #..[i,j,l,m] = aug.k[j,l,i,m], <--> aug.k[i,j,l,m] = ..[l,i,j,m]
    aperm(aug.k, c(3, 1, 2, 4))

  res.arr <- arr.prod +
    #..[i,j,l,m] = arr.prod[j,i, m,l] = h[j,l] * k[i,m]
    aperm(arr.prod, c(2, 1, 4, 3)) -
    #..[i,j,l,m] = arr.prod[i,j,m,l] = h[i,l] * k[j,m]
    aperm(arr.prod, c(1, 2, 4, 3)) -
    #..[i,j,l,m] = arr.prod[j,i,m,l] = h[j,m] * k[i,l]
    aperm(arr.prod, c(2, 1, 3, 4))

  return(res.arr)

}
