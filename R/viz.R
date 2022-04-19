# functions for visualization
# suggested: tidyverse

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
      rows = ggplot2::vars('plt.y'), cols = ggplot2::vars('plt.x'),
      labeller = ggplot2::label_both
    )

  return(plt)

}
