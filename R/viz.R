# functions for visualization
# suggested: tidyverse

approxfunMetric <- function(coord, metric){
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
        'y ~ locfit::lp(%s)', paste(names(dat.coord), collapse = ',')
      )),
      data = w.dat
    )))
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

getMetricEllipsoid <- function(coord, metric, ...){
  # create data.frame of ellipsoid for plotting
  # TBD
  car::ellipse(coord, mat.met, radius = 100, draw = F)
  return(NULL)
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
    )

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


