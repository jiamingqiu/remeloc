# functions for visualization
# suggested: tidyverse

approxfunHD <- function(coord, tensor, bw, method = 'smooth') {
  method <- match.arg(method, c('smooth', 'interp'))
  if(method == 'smooth') {
    res <- approxfunSmooth(coord, tensor, bw)
  } else {
    res <- approxfun_interp2d(coord, tensor)
  }
  return(res)
}

approxfunSmooth <- function(coord, tensor, bw){

  stopifnot(nrow(coord) > 2)

  if(inherits(tensor, 'function'))
    ls.tensor <- apply(coord, 1, tensor, simplify = F)
  else
    ls.tensor <- tensor

  stopifnot(inherits(ls.tensor, 'list'))
  stopifnot(length(ls.tensor) == nrow(coord))
  stopifnot(all(
    sapply(ls.tensor, function(x) all(dim(x) == dim(ls.tensor[[1]])))
  ))

  dim.tensor <- dim(ls.tensor[[1]])
  dd.tensor <- length(ls.tensor[[1]])
  # 1row = 1pnt in coord
  mat.org <- t(sapply(ls.tensor, function(x) as.numeric(x)))

  one.row.f <- function(x) {

    # browser()
    # weighted average to approx inverse, new.emb: new embedded coord

    stopifnot(length(x) == ncol(coord))

    # pseudo Euclidean dist in the embedded space
    pseudo.dist <- colSums( (t(coord) - x) ^ 2 )
    # gaussian weighting
    wt <- dnorm(pseudo.dist, sd = bw)
    wt <- wt / sum(wt)

    use.org <- mat.org[wt > .Machine$double.eps ^ (1/3), , drop = F]
    use.wt <- wt[wt > .Machine$double.eps ^ (1/3)]
    if(length(use.wt) == 0) {
      return(array(0, dim = dim.tensor))
    } else {
      return(
        array(colSums(use.org * use.wt), dim = dim.tensor)
      )
    }

  }

  res.f <- function(target) {
    if(is.data.frame(target))
      x <- as.matrix(target)
    if(!is.matrix(target))
      target <- matrix(target, nrow = 1)

    if(nrow(target) == 1)
      one.row.f(as.numeric(target))
    else
      apply(target, 1, one.row.f, simplify = F)
  }

  return(res.f)

}

approxfun_interp2d <- function(coord, tensor) {
  grid.x <- unique(coord[, 1])
  grid.y <- unique(coord[, 2])
  stopifnot(all(diff(grid.x) > 0) & all(diff(grid.y) > 0))
  stopifnot(identical(
    `colnames<-`(as.matrix(expand.grid(grid.x, grid.y)), NULL),
    coord
  ))
  if(inherits(tensor, 'function')) {
    return(
      approxFun2d(coord, apply(coord, 1, tensor, simplify = F))
    )
  }

  ls.tensor <- tensor
  stopifnot(is.list(ls.tensor))
  dim.tensor <- dim(ls.tensor[[1]])
  dd.tensor <- length(ls.tensor[[1]])
  ls.mat <- lapply(seq(dd.tensor), function(which.element) {
    arr.val <- sapply(ls.tensor, function(tsr) tsr[which.element])
    matrix(arr.val, nrow = length(grid.x))
  })
  one.row.f <- function(x, ...) {
    # browser();QWER
    arr.val <- sapply(seq(dd.tensor), function(which.element) {
      pracma::interp2(
        # remark: pracma seems to iterate y(2nd), unlike R iterate x(1st)
        # see ?pracma::interp2.
        grid.y, grid.x, ls.mat[[which.element]],
        xp = x[2], yp = x[1], ...
      )
    })
    return(array(arr.val, dim = dim.tensor))
  }

  res.f <- function(target, ...) {
    if(is.data.frame(target))
      x <- as.matrix(target)
    if(!is.matrix(target))
      x <- matrix(target, nrow = 1)
    apply(target, 1, function(x, ...) one.row.f(x, ...), simplify = F)
  }

  return(res.f)
}

approxfunTensor <- function(
    coord, tensor, sym.eq = NULL, bw = NULL, optns.locfit = list()
){
  # create approximating (interpolation) function of metric
  # use locfit for now, not good (too smooth), high-d interpolation TBD

  stopifnot(nrow(coord) > 2)

  if(inherits(tensor, 'function'))
    return(
      approxfunTensor(coord, apply(coord, 1, tensor, simplify = F), sym.eq)
    )

  stopifnot(inherits(tensor, 'list'))
  stopifnot(length(tensor) == nrow(coord))
  stopifnot(all(
    sapply(tensor, function(x) all(dim(x) == dim(tensor[[1]])))
  ))

  # dim and rank
  d <- ncol(coord)
  k <- length(dim(tensor[[1]]))

  df.idx <- vecTensorIdx(
    d = d, k = k,
    sym.eq = sym.eq, drop.zero = T
  )
  vec.tsr.basis <- as.matrix(
    df.idx[, stringr::str_detect(names(df.idx), 'tsr\\.e\\d+'), drop = F]
  )


  # compute projection score of tensor to tensor.basis
  # proj.tsr.basis <- tcrossprod(vec.tsr.basis)
  vec.tsr <- sapply(tensor, function(tsr) tsr[df.idx$idx]) # 1 col = 1 tensor
  # sum( (vec.tsr - crossprod(proj.tsr.basis, vec.tsr))^2 )
  vec.tsr <- crossprod(vec.tsr.basis, vec.tsr) #1col = 1 proj.score

  dat.coord <- as.data.frame(coord)
  names(dat.coord) <- sprintf('x%s', seq(d))

  # fit local polynomial
  ls.fit <- apply(vec.tsr, 1, function(arr.comp) {
    w.dat <- dat.coord
    w.dat$y <- arr.comp
    fit <- do.call(locfit::locfit, c(list(
      formula = as.formula(sprintf(
        'y ~ locfit::lp(%s, h = %s)',
        paste(names(dat.coord), collapse = ','),
        deparse(bw)
      )),
      data = w.dat
    ), optns.locfit))
    # fit$call <- NULL
    return(fit)
  })

  # construct return function
  res.f <- function(target){

    # interpolated metric, target: 1point, returns a matrix

    if(!is.matrix(target)) target <- matrix(target, nrow = 1)
    stopifnot(ncol(target) == d)
    # # for the moment
    # stopifnot(nrow(target) > 1)

    new.dat <- as.data.frame(target)
    names(new.dat) <- sprintf('x%s', seq(d))
    n.target <- nrow(new.dat)

    mat.pred <- sapply(ls.fit, function(fit){
      # arr.pred <- rep(NA, n.target)
      arr.pred <- locfit:::predict.locfit(fit, newdata = new.dat)
      stopifnot(length(arr.pred) == n.target)
      return(arr.pred)
    }) # 1col = 1 component func (<~> proj.score), 1row = 1target pnt
    if(n.target == 1) mat.pred <- matrix(mat.pred, nrow = 1)

    # browser();QWER
    mat.tsr <- tcrossprod(vec.tsr.basis, mat.pred) # 1col = 1 target pnt
    ls.res <- apply(mat.tsr, 2, function(arr){
      ps.tsr <- array(0, dim = rep(d, k))
      ps.tsr[df.idx$idx] <- arr
      ps.tsr
    }, simplify = F)

    if(n.target == 1)
      return(ls.res[[1]])
    else
      return(ls.res)

  }

  return(res.f)
}

approxfunMetric <- function(coord, metric, bw = NULL, optns.locfit = list()){
  # create approximating (interpolation) function of metric
  # use locfit for now, not good (too smooth), high-d interpolation TBD

  stopifnot(nrow(coord) > 2)

  # plt.grid <- manifold$genGrid(64)
  if(is.list(metric)){
    ls.metric <- metric
  }else{
    ls.metric <- apply(coord, 1, metric, simplify = F)
  }
  d <- ncol(coord)
  vec.met <- sapply(ls.metric, symMat2Vec) # 1col = 1metric
  dat.coord <- as.data.frame(coord)
  names(dat.coord) <- sprintf('x%s', seq(d))

  # fit local polynomial
  ls.fit <- apply(vec.met, 1, function(arr.comp) {
    w.dat <- dat.coord
    w.dat$y <- arr.comp
    do.call(locfit::locfit, c(list(
      formula = as.formula(sprintf(
        'y ~ locfit::lp(%s, h = %s)',
        paste(names(dat.coord), collapse = ','),
        deparse(bw)
      )),
      data = w.dat
    ), optns.locfit))
  })

  # construct return function
  res.f <- function(target){

    # interpolated metric, target: 1point, returns a matrix

    if(!is.matrix(target)) target <- matrix(target, nrow = 1)
    stopifnot(ncol(target) == d)
    # # for the moment
    # stopifnot(nrow(target) > 1)

    new.dat <- as.data.frame(target)
    names(new.dat) <- sprintf('x%s', seq(d))
    n.target <- nrow(new.dat)

    mat.pred <- sapply(ls.fit, function(fit){
      # arr.pred <- rep(NA, n.target)
      arr.pred <- locfit:::predict.locfit(fit, newdata = new.dat)
      stopifnot(length(arr.pred) == n.target)
      return(arr.pred)
    }) # 1col = 1 component func, 1row = 1target pnt
    if(n.target == 1) mat.pred <- matrix(mat.pred, nrow = 1)

    ls.res <- apply(mat.pred, 1, vec2SymMat, simplify = F)

    if(n.target == 1)
      return(ls.res[[1]])
    else
      return(ls.res)

  }

  return(res.f)
}

#' Draw ellipses of equal geodesic distance to center points
#'
#' @param coord matrix of centers, 1 row is 1 point
#' @param metric list of metric matrices or a function
#' @param radius radius for \code{car::ellipse}, missing then 1
#' @param tissot to provide Tissot indicatrices (equal geodesic distance) or
#' geodesic distance of equal displacement in coordinate chart
#' @param ...
#'
#' @return a data frame whose first column is row number of \code{coord}, while
#' the rest are coordinates.
#'
#' @export
#'
#' @details
#' For a given center point \eqn{p} and metric tensor \eqn{G},
#' if \code{tissot = TRUE}, compute equal geodesic ellipses, points \eqn{x}
#' satisfying \eqn{(x - p)^T G (x - p) = c} with \code{c = radius ^ 2}.
#' Note that this function does not solve for geodesic curves, and the "equal
#' geodesic distance" is computed by local Euclidean approximation, which means
#' one shall not use this for large domain characterization.
#' \cr
#' If \code{tissot = FALSE}, give coordinates for
#' \eqn{(x - p)^T G^{-1} (x - p) = c} with \code{c = radius ^ 2}, which is the
#' cost (geodesic distance) of displacement in coordinate chart. For example,
#' for radius 1, the length of its longer axis equals to the square-root of the
#' larger eigenvalue of the metric tensor, whose square is the geodesic distance
#' of 1 unit displacement in coordinates along such direction.
#'
#' @examples
#' df.ell <- getMetricEllipsoid(rep(0, 2), diag(c(1, 5)), radius = 1)
#' head(df.ell)
#' ggplot2::ggplot(df.ell, ggplot2::aes(x, y, group = idx.pnt)) +
#'   ggplot2::geom_path() + ggplot2::coord_fixed()
#' df.ell <- getMetricEllipsoid(
#'   matrix(seq(4), 2, 2),
#'   list(diag(c(1, 5)), diag(c(1, 1/5))),
#'   radius = 0.5
#' )
#' head(df.ell)
#' ggplot2::ggplot(df.ell, ggplot2::aes(x, y, group = idx.pnt)) +
#'   ggplot2::geom_path() + ggplot2::coord_fixed()
#' # Poincare half-space model for hyperbolic space.
#' df.tissot <- with(spaceHyperbolic(d = 2, model = 'half'), {
#'   getMetricEllipsoid(genGrid(3) + 0.5, metric, radius = 0.1)
#' })
#' df.ell <- with(spaceHyperbolic(d = 2, model = 'half'), {
#'   getMetricEllipsoid(
#'     genGrid(3) + 0.5, metric, radius = 0.1, tissot = F
#'   )
#' })
#' head(df.ell)
#' ggplot2::ggplot(df.tissot, ggplot2::aes(x, y, group = idx.pnt)) +
#'   ggplot2::geom_path() + ggplot2::coord_fixed()
#' ggplot2::ggplot(df.ell, ggplot2::aes(x, y, group = idx.pnt)) +
#'   ggplot2::geom_path() + ggplot2::coord_fixed()
getMetricEllipsoid <- function(coord, metric, radius, tissot = T, ...){
  # create data.frame of ellipsoid for plotting, currently 2-dim only
  # this function gives you ellipses of equal geodesic distance to the coord

  # stopifnot(is.matrix(coord))

  if( !is.data.frame(coord) & !is.matrix(coord) & !is.list(metric) ) {
    coord <- matrix(coord, nrow = 1)
    if(inherits(metric, 'function'))
      metric <- list(metric(coord[1, ]))
    else
      metric <- list(metric)
  }
  if(inherits(metric, 'function'))
    metric <- apply(coord, 1, metric, simplify = F)
  stopifnot(nrow(coord) == length(metric))
  stopifnot(is.list(metric))
  stopifnot(all(sapply(metric, dim) == ncol(coord)))
  if(tissot){
    metric <- lapply(metric, function(x) {
      MASS::ginv(x + 0.01 * max(abs(x)) * diag(ncol(x)))
    })
  }

  if(ncol(coord) != 2) stop('currently only supports d = 2.')

  # browser();QWER
  f.call <- match.call()
  input.args <- as.list(f.call)[-1]
  if(is.null(input.args$radius)) input.args$radius <- 1

  not.my.args <- input.args[setdiff(names(input.args), c(
    'coord', 'metric', 'tissot'
  ))]

  ls.ell <- base::mapply(function(center, mat) {
    # mat <- mat / sum(mat ^ 2) # unify size first
    # mat <- mat
    ls.args <- c(list(
      center = as.numeric(center), shape = mat
    ), not.my.args)
    tm <- as.data.frame(do.call(ellipse2D, ls.args))
    return(tm)
  }, center = asplit(coord, 1), mat = metric, SIMPLIFY = F)

  names(ls.ell) <- seq_along(ls.ell)
  df.ell <- purrr::map_df(ls.ell, ~ .x, .id = 'idx.pnt')
  if(!is.null(colnames(coord))) names(df.ell) <- c('idx.pnt', colnames(coord))
  df.ell$idx.pnt <- as.integer(df.ell$idx.pnt)

  return(df.ell)
}

#' plot graph (edges) on coordinates
#'
#' @param coord n-by-d matrix of coordinates, one row for one point
#' @param edge n-by-2 matrix of edges, one row for one edge
#' @param get.plt.dat \code{FALSE}(default) or \code{TRUE}, to get the plot or
#' data used for plotting.
#' @param max.n.edge max number of edges to show, only when returning plots
#'
#' @return a ggplot2 plot stratified (facet_grid to d-by-d subpanels) by each
#' plotting margin, or a list of data.frames of plotting coordinates and edges.
#'
#' @details The returned list (when \code{get.plt.dat == TRUE}) contains two
#' elements: \code{coord} and \code{edge}, each containing the coordinates
#' \code{x, y, xend, yend} used to draw the plot. Also, \code{idx.pnt} and
#' \code{idx.edge} columns record the row number of the plotting point or edge
#' in the input \code{coord} or \code{edge}.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' d <- 3
#' coord <- matrix(runif(d * 50), ncol = d)
#' edge <- allEdge(coord, 0.25)
#' plotGraph(coord, edge)
plotGraph <- function(coord, edge, get.plt.dat = F, max.n.edge = 500){

  # use ggplot to plot a graph (defined by edge) on coord

  # DEV ONLY
  # coord <- matrix(0, nrow = 5, ncol = 2)
  # coord[-1, ] <- as.matrix(expand.grid(v1 = c(-1, 1), v2 = c(-1, 1)))
  # edge <- allEdge(coord, 1)
  # set.seed(1)
  # d <- 3
  # coord <- matrix(runif(d * 50), ncol = d)
  # edge <- allEdge(coord, 0.25)
  # END DEV ONLY

  idx.plt.margin <- utils::combn(ncol(coord), 2)
  idx.plt.margin <- apply(idx.plt.margin, 2, function(x) x, simplify = F)
  plt.coord <- lapply(idx.plt.margin, function(plt.margin){
    # subset the plotting margin (column)
    coord[, plt.margin]
  })
  plt.edge <- lapply(plt.coord, function(x) {
    tm.mat <- matrix(0, nrow = nrow(edge), ncol = 2 * ncol(x))
    tm.mat[, seq(ncol(x))] <- x[edge[, 1], ]
    tm.mat[, seq(ncol(x)) + ncol(x)] <- x[edge[, 2], ]
    # tm.mat col: x, y, xend, yend
    return(tm.mat)
  })

  plt.edge <- purrr::pmap(list(idx.plt.margin, plt.edge), ~ {
    # browser();QWER
    df <- as.data.frame(..2)
    names(df) <- c('x', 'y', 'xend', 'yend')
    df$idx.edge <- seq(nrow(df))
    df$plt.x <- ..1[1]
    df$plt.y <- ..1[2]
    df$plt.margin <- paste(..1, collapse = '-')
    return(df)
  })
  plt.edge <- dplyr::bind_rows(plt.edge)
  plt.coord <- purrr::pmap(list(idx.plt.margin, plt.coord), ~ {
    df <- as.data.frame(..2)
    names(df) <- c('x', 'y')
    df$idx.pnt <- seq(nrow(df))
    df$plt.x <- ..1[1]
    df$plt.y <- ..1[2]
    df$plt.margin <- paste(..1, collapse = '-')
    return(df)
  })
  plt.coord <- dplyr::bind_rows(plt.coord)

  if(get.plt.dat){
    return(list(
      coord = plt.coord, edge = plt.edge
    ))
  }

  if(nrow(edge) > max.n.edge){
    # thinning
    warning(sprintf(
      'total edges: %s, only showing %s.', nrow(edge), max.n.edge
    ))
    idx.keep <- sample.int(nrow(edge), size = max.n.edge)
    plt.edge <- split(plt.edge, plt.edge$plt.margin)
    plt.edge <- lapply(plt.edge, function(df) {
      df <- df[idx.keep, ]
      # df$idx.edge <- seq(nrow(df))
      return(df)
    })
    plt.edge <- dplyr::bind_rows(plt.edge)

    plt.coord <- plt.edge %>% dplyr::select(-c('x', 'y')) %>%
      dplyr::rename('x' = 'xend', 'y' = 'yend')
    plt.coord <- dplyr::bind_rows(
      plt.coord, plt.edge %>% dplyr::select(-c('xend', 'yend'))
    )
    plt.coord <- dplyr::distinct(plt.coord)
  }

  plt <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = plt.edge, mapping = ggplot2::aes(
      x = x, y = y, xend = xend, yend = yend
      # , color = as.integer(as.factor(idx.edge))
    ), show.legend = F) +
    ggplot2::geom_point(
      data = plt.coord,
      mapping = ggplot2::aes(x = x, y = y),
      color = 'red', inherit.aes = FALSE
    ) +
    # ggplot2::scale_color_gradientn(colours = hcl.colors(10, 'Zissou 1')) +
    ggplot2::facet_grid(
      # rows = ggplot2::vars('plt.y'), cols = ggplot2::vars('plt.x'),
      as.formula('plt.y ~ plt.x'),
      labeller = ggplot2::label_both
    ) #+ ggplot2::coord_fixed() # save for user, in case sf?

  return(plt)

}


#' get data frame from metric for plotting
#'
#' @param coord a matrix of targeted points
#' @param metric metric function or a list of metric at \code{coord}
#'
#' @return a data.frame
#' @export
#'
#' @examples
#' df <- with(spaceHyperbolic(d = 2, model = 'half'), {
#'   getMetricDF(genGrid(10) + 1, metric)
#' })
#' ggplot2::ggplot(df) +
#'   ggplot2::geom_tile(ggplot2::aes(x = x1, y = x2, fill = value)) +
#'   ggplot2::facet_wrap(~ component, labeller = ggplot2::label_both) +
#'   ggplot2::coord_fixed()
getMetricDF <- function(coord, metric){

  # get a data.frame for metric at coord, metric can be a list of matrices

  # plt.grid <- manifold$genGrid(64)
  if(is.list(metric)){
    ls.metric <- metric
  }else{
    ls.metric <- apply(coord, 1, metric, simplify = F)
  }
  d <- ncol(coord)

  stopifnot(all(
    sapply(ls.metric, function(x) identical(dim(x), rep(d, 2)) )
  ))

  df.metric <- t(sapply(ls.metric, symMat2Vec)) # one row = one point
  df.metric <- as.data.frame(df.metric)
  # add coordinates
  df.metric <- cbind(coord, df.metric)
  # names

  idx.vec <- vecQuadIdx(d)
  names(df.metric) <- c(
    sprintf('x%s', seq(d)),
    sprintf('dx%sdx%s', idx.vec[1, ], idx.vec[2, ])
  )

  df.metric <- tidyr::pivot_longer(df.metric, cols = tidyselect::all_of(
    sprintf('dx%sdx%s', idx.vec[1, ], idx.vec[2, ])
  ), names_to = 'component', values_to = 'value')

  return(df.metric)

}


ellipse2D <- function (center, shape, radius, segments = 51) {

  # modified based on car::ellipse, give coordinates for ellipse
  # {x: t(x - center) %*% solve(shape) %*% (x - center) = radius ^ 2}

  if (!(is.vector(center) && 2 == length(center)))
    stop("center must be a vector of length 2")
  if (!(is.matrix(shape) && all(2 == dim(shape))))
    stop("shape must be a 2 by 2 matrix")
  if (max(abs(shape - t(shape)))/max(abs(shape)) > 1e-10)
    stop("shape must be a symmetric matrix")

  angles <- (0:segments) * 2 * pi/segments
  unit.circle <- cbind(cos(angles), sin(angles))
  Q <- chol(shape, pivot = TRUE)
  order <- order(attr(Q, "pivot"))
  ellipse <- t(center + radius * t(unit.circle %*% Q[, order]))

  colnames(ellipse) <- c("x", "y")
  return(ellipse)
}



