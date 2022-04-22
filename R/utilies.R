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
#' @examples
#' set.seed(1)
#' d <- 3
#' coord <- matrix(runif(d * 50), ncol = d)
#' graph <- genGraph(10, coord, local.reach = 0.5)
#' plotGraph(coord, graph)
genGraph <- function(n, coord, local.reach){

  if(!is.matrix(coord))
    stop('need more than 1 points.')
  stopifnot(all(local.reach >= 0))
  stopifnot(length(local.reach) == 1 | length(local.reach) == ncol(coord))
  stopifnot(n > 0)

  coord.pnts <- coord

  # index of all possible pairs within range of local.reach
  all.edge <- allEdge(coord, local.reach)
  all.edge <- all.edge[all.edge[, 1] != all.edge[, 2], ]
  res <- all.edge[
    sample(nrow(all.edge), size = n), , drop = F
  ]

  return(res)
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
#' Note that self-loop will be included.
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

# the R version is faster when there is larger number of points
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

#' change format of a graph (edge)
#'
#' @param graph a list or matrix of edges of the graph
#' @param format \code{"matrix"} or \code{"list"}, the target format, if
#' missing, will toggle the input type, i.e., change \code{matrix} input to
#' \code{list}, and vise versa
#' @param ... additional argument passed to \code{edgeMat2List}
#'
#' @return the reformatted edges of the graph
#' @export
#' @family {loc_graph}
#' @examples
#' set.seed(1)
#' d <- 2
#' coord <- matrix(runif(d * 10), ncol = d)
#' edge.mat <- locWindow(c(0, 0), coord, 0.5, 'matrix')
#' edge.list <- locWindow(c(0, 0), coord, 0.5, 'list')
#' identical(formatGraph(edge.mat), edge.list)
#' identical(formatGraph(edge.list), edge.mat)
formatGraph <- function(graph, format, ...){

  stopifnot( inherits(graph, c('matrix', 'list')) )
  if(missing(format)){
    format <- setdiff(c('matrix', 'list'), class(graph))
  }else{
    format <- match.arg(format, c('matrix', 'list'))
    if(inherits(graph, format)) return(graph)
  }

  if(format == 'matrix')
    return(edgeList2Mat(graph))
  else
    return(edgeMat2List(graph, ...))
}

#' translate a matrix of graph edges into list
#'
#' @param edge.mat a 2-col matrix of graph edges (pair of indices)
#' @param names_from which of the 2 column to go to list name,
#' either integer index or character for column name.
#'
#' @return a named list whose elements are indices of nodes connected to the
#' node corresponding to the list name.
#'
#' @export
#' @family {loc_graph}
#' @examples
#' set.seed(1)
#' d <- 2
#' coord <- matrix(runif(d * 10), ncol = d)
#' edge.mat <- locWindow(c(0, 0), coord, 0.5, 'matrix')
#' edge.list <- locWindow(c(0, 0), coord, 0.5, 'list')
#' identical(edgeMat2List(edge.mat), edge.list)
edgeMat2List <- function(edge.mat, names_from = NULL){
  # translate a matrix of graph edges into list
  # basically a wrapper around split and build-your-own-wheel tidyr::nest

  stopifnot(is.matrix(edge.mat))
  stopifnot(ncol(edge.mat) == 2)

  # if( length(unique(edge.mat[, 1])) > length(unique(edge.mat[, 2])) ) {
  #   # if the 2nd col is the one with fewer levels
  #   edge.mat <- edge.mat[, c(2, 1)]
  # }

  if(is.null(colnames(edge.mat)))
    colnames(edge.mat) <- seq(2)
  if(is.null(names_from)) names_from <- colnames(edge.mat)[1]
  idx.name <- which(names_from == colnames(edge.mat))
  # browser();QWER
  tm <- edge.mat[, -idx.name]
  names(tm) <- NULL
  return(split(tm, edge.mat[, idx.name]))

}

#' translate a matrix of graph edges into list
#'
#' @param edge.list a named list for edges (indices)
#'
#' @return a 2-col matrix of edges (pair of indices)
#' @export
#'
#' @family {loc_graph}
#'
#' @examples
#' set.seed(1)
#' d <- 2
#' coord <- matrix(runif(d * 10), ncol = d)
#' edge.mat <- locWindow(c(0, 0), coord, 0.5, 'matrix')
#' edge.list <- locWindow(c(0, 0), coord, 0.5, 'list')
#' identical(edgeList2Mat(edge.list), edge.mat)
edgeList2Mat <- function(edge.list){
  # translate a list of graph edges into edge

  if(is.null(names(edge.list))) {
    warning('unnamed list input, using seq_along(edge.list).')
    names(edge.list) <- seq_along(edge.list)
  }
  # make into a matrix
  tm <- unlist(edge.list)
  mat <- matrix(0L, nrow = length(tm), ncol = 2) # 0L to keep integer type
  mat[, 1] <-
    rep(as.integer(names(edge.list)), times = sapply(edge.list, length))
  mat[, 2] <- tm

  return(mat)

}

#' compute difference in coordinates corresponding to local graph
#'
#' @param coord a matrix of coordinates (1 row = 1 point)
#' @param edge a matrix or list defining the edges connecting \code{coord}
#' @param return.format the format for return, \code{"keep"} or
#' \code{"full.matrix"}
#'
#' @return a matrix or list of difference in coordinates
#'
#' @details
#' If \code{return.format == 'keep'}, will follow input format of \code{edge},
#' say, if \code{edge} is a matrix, then return a matrix of difference in the
#' coordinates of edges (1st column - 2nd column in \code{edge}), with one row
#' corresponds to one edge.
#' Otherwise, \code{edge} is a list, then return a list of such matrices.
#' \cr
#' If \code{return.format == 'full.matrix'}, return a matrix of
#' \code{2 + ncol(coord)} columns, where the first 2 are matrix of \code{edge},
#' while the last few are difference in coordinates.
#'
#' @export
#'
#' @examples TBD
deltaCoord <- function(coord, edge, return.format = 'keep'){
  # compute difference in coordinates according to edge
  edge.in <- edge
  if(typeof(edge.in) == 'list') {
    edge <- edgeList2Mat(edge)
  }
  stopifnot(ncol(edge) == 2)
  stopifnot(is.matrix(coord))
  stopifnot(min(edge) >= 1 & max(edge) <= nrow(coord))

  d <- ncol(coord)
  coord.diff <- coord[edge[, 1], , drop = F] - coord[edge[, 2], , drop = F]

  return.format <- match.arg(return.format, c('keep', 'full.matrix'))

  if(return.format == 'keep'){
    if(typeof(edge.in) != 'list') {
      return(coord.diff)
    }else{
      coord.diff <- split(coord.diff, edge[, 1])
      coord.diff <- lapply(coord.diff, function(x) matrix(x, ncol = d))
      return(coord.diff)
    }
  }else{
    res <- matrix(0, nrow = nrow(edge), ncol = 2 + d)
    res[, seq(2)] <- edge
    res[, 2 + seq(d)] <- coord.diff
    colnames(res) <- c('p1', 'p2', sprintf('coordd.%s', seq(d)))
    return(res)
  }

  stop('something wrong.')

}

#' identify local neighbors of target points
#'
#' @param target a matrix for coordinates of targets, 1 row = 1 point
#' @param coord a matrix for coordinates of points, 1 row = 1 point
#' @param local.reach approximate size of local neighborhood square
#' (treated as Euclidean space), recycled if necessary
#' @param return \code{"list"} (default) / \code{"matrix"}, return format
#'
#' @return a list or matrix of indices of \code{coord} that is within
#' neighborhood of \code{target}.
#'
#' @details If \code{return == "list"}, then a list with length equal to
#' \code{nrow(target)}, whose ith elements are indices (row) for points in
#' \code{coord} that are within \code{local.reach} of \code{target[i, ]}.
#' \cr
#' If \code{return == "matrix"} then a 2-column matrix where first column are
#' indices of points in \code{target} and second are indices of \code{coord}
#' that is within \code{local.reach} of the target point.
#' \cr
#' Note that this function only consider a square neighborhood of coordinate
#' values, no geometry is involved.
#'
#' @export
#'
#' @seealso edgeMat2List for changing format
#'
#' @examples set.seed(1)
#' d <- 2
#' coord <- matrix(runif(d * 10), ncol = d)
#' locWindow(c(0, 0), coord, 0.5, 'matrix')
locWindow <- function(target, coord, local.reach, return = 'list'){
  # a function finding neighbors
  # args:
  #   coord: a matrix for coordinates of points, 1 row = 1 point
  #   width: width of intervals
  #   return: list(default)/data.frame, what to return
  # returns: a list/data.frame of 2 column with elements being seq(nrow(coord)),
  #   for example, a row of (1, 2) means the point coord[2, ] is within
  #   the box of width centering at point coord[1, ].
  # Details: only judging by difference in coordinates, does not care about
  #   geometry. Also, the ith list will start from i.

  # target <- matrix(rnorm(2 * 10), ncol = 2, nrow = 10)
  # coord <- matrix(rnorm(2 * 5000), ncol = 2, nrow = 5000)
  # local.reach <- 0.15

  stopifnot(
    all(is.finite(coord)) & all(is.finite(local.reach)) & all(is.finite(target))
  )

  if(!is.matrix(coord))
    coord <- matrix(coord, nrow = 1)
  if(!is.matrix(target))
    target <- matrix(target, nrow = 1)
  d <- ncol(coord)
  n <- nrow(target)

  stopifnot(all(local.reach >= 0))
  stopifnot(length(local.reach) == 1 | length(local.reach) == d)
  stopifnot(ncol(target) == d)
  if(length(local.reach) == 1)
    local.reach <- rep(local.reach, d)

  t.coord <- t(coord) # one col for one point, for vectorization
  tm.res <- lapply(seq(n), function(idx){
    which(
      colSums(
        abs(t.coord - target[idx, ]) <= local.reach
      ) == d
    )
  })
  names(tm.res) <- seq_along(tm.res)

  return <- match.arg(return, c('list', 'matrix'))
  if(return == 'list'){
    return(tm.res)
  }
  tm <- unlist(tm.res)
  mat <- matrix(0L, nrow = length(tm), ncol = 2) # 0L to keep integer type
  # note that if length(element) is 0, rep times = 0, so we are fine.
  mat[, 1] <- rep(seq_along(tm.res), times = sapply(tm.res, length))
  mat[, 2] <- tm
  return(mat)
}
