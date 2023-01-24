# utilities

inConvexPolygon <- function(p, vertex) {
  # determine if p is inside the convex hull of vertex
  # p & vertex: 1 row for 1 point
  # Math unvalidated yet ,use with caution

  stopifnot(is.matrix(vertex))
  # stopifnot(length(p) == ncol(vertex))

  hull <- geometry::convhulln(vertex)
  return(
    geometry::inhulln(hull, p)
  )

}

# do.call.tommy <- function(what, args, ...) {
#   # credit to Tommy at
#   # https://stackoverflow.com/questions/10022436/do-call-in-combination-with
#   if(is.character(what)){
#     fn <- strsplit(what, "::")[[1]]
#     what <- if(length(fn)==1) {
#       get(fn[[1]], envir=parent.frame(), mode="function")
#     } else {
#       get(fn[[2]], envir=asNamespace(fn[[1]]), mode="function")
#     }
#   }
#
#   do.call(what, as.list(args), ...)
# }
# do.call.jeroen <- function(what, args, ...){
#
#   # credit to Jeroen at
#   # https://stackoverflow.com/questions/10022436/do-call-in-combination-with
#
#   if(is.function(what)){
#     what <- deparse(as.list(match.call())$what);
#   }
#   myfuncall <- parse(text=what)[[1]];
#   mycall <- as.call(c(list(myfuncall), args));
#   eval(mycall, ...);
# }

relFrob <- function(esti, true, eps = .Machine$double.eps^(1/3)){
  # relative Frobenius error
  frob.true <- sqrt(sum(true ^ 2))
  frob.err <- sqrt(sum( (esti - true) ^ 2 ))
  if(frob.true < eps)
    return(frob.err)
  else
    return(frob.err / frob.true)
}

# combination number encode/decoder
comb2_vec2idx <- function(vec, n){
  # vec-th combination to the actual combination (i, j)
  lookup <- c(0L, cumsum(rev(seq(n - 1))))
  res <- matrix(0L, ncol = length(vec), nrow = 2)
  res[1, ] <- sapply(vec, function(x) max(which(x > lookup)) )
  res[2, ] <- vec - lookup[res[1, ]]
  res[2, ] <- (seq(n))[res[2, ] + res[1, ]]
  return(res)
}
comb2_idx2vec <- function(idx, n){
  # order of combination (i, j) in combn(n, 2) without invoking combn
  # (n - 1) + (n - 2) + ... + (n - i + 1) + j - i
  if(!is.matrix(idx)) idx <- as.matrix(idx)
  stopifnot(nrow(idx) == 2)
  stopifnot(min(idx) == 1 & max(idx) == n)
  return(
    ((2L * n - idx[1, ]) * (idx[1, ] - 1L)) / 2L + idx[2, ] - idx[1, ]
  )
}


#' Reindex array
#'
#' @description Rearrange array so that \code{x_out[to] == x_in[from]}.
#' Just a wrapper around \code{aperm} for easier use.
#'
#' @param x an array
#' @param to new subscripts
#' @param from old subscripts
#'
#' @return permuted \code{x}
#' @export
#'
#' @examples
#' reindex(matrix(seq(4), 2, 2), 'ij', 'ji')
#' tst.arr <- array(seq(16), rep(2, 4))
#' reindex(tst.arr, 'ijkl', 'jilk')
#' for(i in seq(2)){ for(j in seq(2)){
#'   for(k in seq(2)){ for(l in seq(2)){
#'     stopifnot(identical(
#'       reindex(tst.arr, 'ijkl', 'jilk')[i, j, k, l],
#'       tst.arr[j, i, l, k]
#'     ))
#' }}}}
reindex <- function(x, to, from){

  # a wrapper around aperm for easier use
  dim.x <- dim(x)
  idx.new <- strsplit(to, '')[[1]]
  idx.old <- strsplit(from, '')[[1]]
  stopifnot(length(idx.new) == length(dim.x) & length(idx.old) == length(dim.x))
  stopifnot(identical(sort(idx.new), sort(idx.old)))

  perm <- match(idx.new, idx.old)

  return(aperm(x, perm))

}

#' Index points
#' @description create a look-up dictionary for input points, predominantly used
#' to simplify data structure of graphs.
#' @param pnts a matrix of coordinates of points
#'
#' @return a list of the look-up dictionary matrix \code{coord} and an integer
#' vector \code{idx} indexing the input \code{pnts}.
#'
#' @export
#'
#' @examples
#' pnts <- matrix(seq(9), ncol = 3)
#' pnts <- rbind(pnts, c(1, 4, 7), matrix(rnorm(3 * 2), ncol = 3), c(2, 5, 8))
#' indexPoints(pnts)
indexPoints <- function(pnts){
  # index points for easier edge indexing and return a list

  # DEV ONLY
  # pnts <- matrix(seq(9), ncol = 3)
  # pnts <-
  #   rbind(pnts, c(1, 4, 7), matrix(rnorm(3 * 100), ncol = 3), c(2, 5, 8))
  # END DEV ONLY

  if(!is.matrix(pnts)) pnts <- matrix(pnts, nrow = 1)
  coord <- pnts
  idx.unique <- !base::duplicated(coord, MARGIN = 1, )
  coord <- coord[idx.unique, , drop = F]
  # if unique, idx + 1, otherwise still
  res.idx <- cumsum(idx.unique)
  # match those not unique
  dup.idx <- match(
    apply(pnts[!idx.unique, , drop = F], 1, paste, collapse = ' '),
    apply(coord, 1, paste, collapse = ' ')
  )

  res.idx[!idx.unique] <- dup.idx
  return(list(
    idx = res.idx, coord = coord
  ))
}
# microbenchmark::microbenchmark(list = alist(
#   string = match( # faster
#     apply(pnts[!idx.unique, , drop = F], 1, paste, collapse = ' '),
#     apply(coord, 1, paste, collapse = ' ')
#   ),
#   asplit = match(
#     asplit(pnts[!idx.unique, , drop = F], 1),
#     asplit(coord, 1)
#   )
# ))


#' Random generator for a local graph
#'
#' @param n number of edge to generate
#' @param coord matrix of coordinates, one row for one point
#' @param local.reach approximate size of local neighborhood
#' (treated as Euclidean space)
#' @param exact whether to call \code{allEdge} and select from all local edges
#'
#' @return a n-by-2 matrix for indices of edge, one row for one edge.
#' @details will based on coordinates, i.e., no edge if coordinates are too far
#' apart, in the sense that difference in coord is greater than \code{range}.
#' @export
#'
#' @examples
#' set.seed(1)
#' d <- 2
#' coord <- matrix(runif(d * 50), ncol = d)
#' graph <- genGraph(10, coord, local.reach = 0.5)
#' plotGraph(coord, graph)
#' plotGraph(coord, genGraph(10, coord, local.reach = 0.5, T))
#' plotGraph(coord, genGraph(100, coord, local.reach = 0.2))
genGraph <- function(n, coord, local.reach, exact = F){

  # browser();QWER
  if(!is.matrix(coord))
    stop('need more than 1 points.')
  stopifnot(all(local.reach >= 0))
  stopifnot(length(local.reach) == 1 | length(local.reach) == ncol(coord))
  stopifnot(n > 0)

  coord.pnts <- coord

  if(exact){
    # index of all possible pairs within range of local.reach
    all.edge <- allEdge(coord, local.reach)
    all.edge <- all.edge[all.edge[, 1] != all.edge[, 2], ]
    if(nrow(all.edge) < n) warning(
      'fewer candidate than requested, try larger local.reach'
    )
    res <- all.edge[
      sort(sample.int(nrow(all.edge), size = min(n, nrow(all.edge)))), ,
      drop = F
    ]

  }else{

    max.iter <- 100
    iter <- 0
    res <- matrix(nrow = 0, ncol = 2)

    while(nrow(res) < n & iter < max.iter){

      iter <- iter + 1

      # store index for later recovery
      use.idx <- sort(sample.int(nrow(coord), size = min(10^4, nrow(coord))))

      use.coord <- coord[use.idx, , drop = F]
      suppressWarnings({
        use.edge <-
          genGraph(n = n, use.coord, local.reach = local.reach, exact = T)
      })

      # translate index
      use.edge[, 1] <- use.idx[use.edge[, 1]]
      use.edge[, 2] <- use.idx[use.edge[, 2]]

      res <- rbind(res, use.edge)

      res <- res[
        !base::duplicated(sprintf('%s - %s', res[, 1], res[, 2])), , drop = F
      ]
      res <- res[order(res[, 1]), ]

    }
    # root.pnt <- apply(coord, 2, function(x){
    #   seq(min(x) + local.reach / 2, max(x), by = local.reach / 2)
    # }, simplify = F)
    #
    # root.pnt <- as.matrix(do.call(expand.grid, root.pnt))
    #
    # # plot(coord)
    # # points(root.pnt, col = 'red')
    #
    # ls.nbhd <- locWindow(root.pnt, coord, local.reach / 2, return = 'list')
    # # drop those with fewer than 2 points
    # ls.nbhd <- ls.nbhd[sapply(ls.nbhd, length) >= 2]
    # # sample proportion in each local window
    # sample.prop <-
    #   n / sum(sapply(ls.nbhd, function(x) length(combn(length(x), 2))))
    # res <- lapply(ls.nbhd, function(x){
    #   idx <- combn(length(x), 2)
    #   idx <- idx[, sort(sample.int( # sort for uniqueness
    #     ncol(idx), size = min(ceiling(ncol(idx) * sample.prop * 3), ncol(idx))
    #   )), drop = F]
    #   return(cbind(x[idx[1, ]], x[idx[2, ]]))
    # })
    # res <- do.call(rbind, res)
    #
    # size
    if(nrow(res) < n) stop(
      'fewer candidate than requested, try exact method or larger local.reach'
    )
    res <- res[sort(sample.int(nrow(res), size = n)), , drop = F]

    # sort
    res <- res[order(res[, 1]), ]
  }

  return(res)
}

#' Random generate for edges comparison in local graph
#'
#' @param n number of comparison
#' @param edge a list or matrix of edges
#'
#' @return a n-by-4 matrix
#' @export
#'
#' @details Uses a two-phase generation scheme, where we pick the primary
#' edges, then the secondary ones are selected by looking at nodes connected to
#' the primary ones.
#' \code{n} can be length 1 or 2, if length 2, the second number is used in
#' second phase for multiple comparison. It is possible that the resulting
#' sample is fewer than requested (especially with sparse graph), in which case
#' one should consider use larger \code{n[2]}.
#' \cr
#' Besides, for unique labeling edges in the resulting matrix, its first column
#' will be smaller than the second column; the third column smaller than the
#' fourth and the first column; and the second column smaller than the fourth
#' column.
#'
#' #' @examples
#' d <- 3
#' manifold <- spaceEuclidean(d)
#' set.seed(1)
#' obsv.coord <- manifold$genPnt(10)
#' all.edge <- allEdge(obsv.coord, local.reach = 0.5)
#' all.edge <- all.edge[all.edge[, 1] < all.edge[, 2], ]
#' plotGraph(obsv.coord, all.edge)
#' genCompareGraph(10, all.edge)
genCompareGraph <- function(n, edge, coord, local.reach){

  if(length(n) == 1)
    n <- c(n, 1)
  stopifnot(all(n >= 1))
  # genenrate graph comparing edges
  mat.edge <- formatGraph(edge, format_to = 'matrix')
  mat.edge <- mat.edge[mat.edge[, 1] != mat.edge[, 2], , drop = F]
  mat.edge <- t(apply(mat.edge, 1, sort)) # for later uniqueness

  if(missing(coord) & missing(local.reach)){

    # if sampling only using edge specification, slow

    tm <- formatGraph(mat.edge, format_to = 'list')
    ls.edge <- rep(list(NULL), max(mat.edge))
    names(ls.edge) <- seq_along(ls.edge)
    ls.edge[names(tm)] <- tm

    primary.edge <- mat.edge[
      sort(sample.int(
        nrow(mat.edge) - 1,
        size = min(nrow(mat.edge) - 1, 2 * n[1])
      )), , drop = F
    ]

    secondary.edge <- lapply(asplit(primary.edge, 1), function(edge.1){
      # browser();QWER
      # # compare randomly others, could be far away
      # edge.option <- ls.edge[seq(edge.1[1], length(ls.edge))]
      # compare nearby
      node.option <- setdiff(unique(unlist(ls.edge[edge.1])), edge.1)
      edge.option <- ls.edge[node.option]
      edge.option <- lapply(edge.option, function(x) x[x > edge.1[2]])
      edge.option <- formatGraph(edge.option, format_to = 'matrix')
      n.option <- nrow(edge.option)
      edge.2 <- edge.option[
        sample.int(n.option, size = min(n.option, n[2])), , drop = F
      ]
      res <- matrix(0L, ncol = 4, nrow = nrow(edge.2))
      res[, seq(2)] <- rep(edge.1, each = nrow(res))
      res[, 2 + seq(2)] <- edge.2
      return(res)
    })

    res <- do.call(rbind, secondary.edge)
    if(nrow(res) < prod(n))
      warning('limited choices, consider increase n[2].')
    res <- res[sort(sample.int(
      nrow(res), size = min(nrow(res), prod(n)))
    ), , drop = F]
    return(res)

  }else{

    # if also providing additional information, can be faster
    # ad hoc coordinates for edges

    pseudo.coord <- (coord[mat.edge[, 1], ] + coord[mat.edge[, 2], ]) / 2

    # too slow
    # pseudo.edge <- allEdge(pseudo.coord, local.reach = local.reach)
    # pseudo.edge <- pseudo.edge[pseudo.edge[, 1] != pseudo.edge[, 2], , drop = F]
    # pseudo.edge <-
    #   pseudo.edge[sample.int(nrow(pseudo.edge), prod(n)), , drop = F]
    pseudo.edge <- genGraph(prod(n), pseudo.coord, local.reach = local.reach)
    res <- cbind(mat.edge[pseudo.edge[, 1], ], mat.edge[pseudo.edge[, 2], ])

    return(res)

  }
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
    tm.res <- allEdge_cpp(
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
#   res <- allEdge_cpp(
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
#' @param format_to \code{"matrix"} or \code{"list"}, the target format, if
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
formatGraph <- function(graph, format_to, ...){

  stopifnot( inherits(graph, c('matrix', 'list')) )
  if(missing(format_to)){
    format_to <- setdiff(c('matrix', 'list'), class(graph))
  }else{
    format_to <- match.arg(format_to, c('matrix', 'list'))
    if(inherits(graph, format_to)) return(graph)
  }

  if(format_to == 'matrix')
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
