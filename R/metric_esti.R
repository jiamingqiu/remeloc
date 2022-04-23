### core functions for metric estimation


##### options setting ##########################################################

# a handling one, another guessing model, and one drop unnecessary ones for glm

idx_glm_optns <- function(optns){
  # return index of options list that is passed to glm.fit
  # by dropping what our method uses.
  res <- which(!(names(optns) %in% c(
    'n.local', 'method.trim', 'metric.only',
    'local.reach', 'intercept', 'type'
  )))
  return(res)
}

#' set or update computational options
#'
#' @param optns a list of options
#' @param ... what to set or update
#'
#' @return a list of options
#' @export
#'
#' @details Possible options are
#' \describe{
#'    \item{\code{method.trim}}{
#'        \code{"proximity"}(default)/\code{"random"}, method to select local
#'        responses, either random, or by proximity to target points. Here
#'        distance is taken assuming points are in Euclidean space regardless of
#'        truth.
#'    }
#'    \item{\code{n.local}}{
#'        max and min number of local responses to use. Missing then default to
#'        \code{c(10, 1000)}. If provided only one value, will use as max, and
#'        min will be set to 1.
#'    }
#'    \item{\code{metric.only}}{
#'        \code{TRUE}(default)/\code{FALSE}, whether to return only the
#'        estimated local metric tensor per target point, or to also include
#'        additional results, such as local observations, fits, and options.
#'    }
#' }
#'
#' @family {locMetric}
#'
#' @examples
#' set_optns()
#' set_optns(method.trim = 'random')
set_optns <- function(optns = list(), ...){
  # setting/updating options
  # browser();QWER

  # fill in the input keys
  f.call <- match.call()
  ls.args <- as.list(f.call)[-1]
  ls.args <- ls.args[names(ls.args) != 'optns']
  optns[names(ls.args)] <- ls.args

  ## fill in the defaults
  # if locally more data available how to trim down the size
  optns$method.trim <- match.arg(optns$method.trim, c(
    'proximity', 'random'
  ))

  # max and min number of local responses to use
  if(is.null(optns$n.local)) optns$n.local <- c(10, 1000)
  stopifnot(all(optns$n.local >= 1))
  if(length(optns$n.local) == 1) optns$n.local <- c(1, optns$n.local)

  # what to return after estimation
  if(is.null(optns$metric.only)) optns$metric.only <- TRUE
  stopifnot(inherits(optns$metric.only, 'logical'))

  return(optns)

}

#' guess model based on input
#'
#' @param optns a list of options
#' @param formula model formula, optional
#' @param data data, optional
#'
#' @return a list with model options, computational options will also be set.
#' Note that guess will only take place when necessary options are not specified
#' in \code{optns}.
#'
#' @export
#'
#' @details Possible options are
#' \describe{
#'   \item{\code{type}}{
#'     \code{"distance"}, \code{"thresholding"}, or \code{"compare"}.
#'     If \code{data} is supplied, will use the number of unique
#'     values in its first column to determine whether it is binary response,
#'     then, the number of its columns used to determine how many edges are
#'     involved.
#'     If only \code{formula} is supplied, will choose from \code{"distance"} or
#'     \code{"compare"}, depending on number of edges deduced from
#'     \code{formula}.
#'     If non of \code{formula} and \code{data} not supplied, default to
#'     \code{"distance"}.
#'   }
#'   \item{\code{intercept}}{
#'     \code{FALSE} or \code{TRUE}, whether to add an intercept term.
#'     Will default to \code{FALSE} unless \code{type = "thresholding"}.
#'   }
#'   \item{\code{family}}{
#'     passed to \code{glm.fit}, Default:
#'     missing for \code{type = "distance"};
#'     \code{binomial(link = "logit")} for
#'        \code{type = c("thresholding", "compare")}.
#'     Note that if no \code{family} specified in \code{type = "distance"}, will
#'     use a simple QR-solve for OLS instead of \code{glm.fit}.
#'   }
#'   \item{\code{...}}{
#'     passed to \code{glm.fit}, say, specify weight when using
#'     \code{glm.fit}.
#'   }
#' }
#'
#'
#' @family {locMetric}
#'
#' @examples
#' guess_model()
guess_model <- function(optns = list(), formula, data, ...){
  # a function used to guess model based on formula or data then fill-in
  # unspecified optns

  # fill-in default computational options first
  optns <- set_optns(optns, ...)

  # if(missing(formula)) stopifnot(!missing(data))
  # if(missing(data)) stopifnot(!missing(formula))

  # then guess
  guess.optns <- list()
  if(missing(formula) & missing(data)){

    guess.optns$type <- 'distance'

  }else if(missing(formula)){

    stopifnot(!missing(data))
    if( length(unique(data[, 1])) > 2 ){
      if((ncol(data) - 1) %% 2 != 0) stop(
        'continuous response, but with edge comparison format.'
      )
      guess.optns$type <- 'distance'
    }else if( (ncol(data) - 1) %% 3 == 0 ){
      guess.optns$type <- 'thresholding'
    }else{
      guess.optns$type <- 'compare'
    }

  }else{

    stopifnot(!missing(formula))
    n.edges <- 1 + length(stringr::str_extract_all(
      deparse(formula, width.cutoff = 200L), ':'
    )[[1]])
    if(n.edges == 2){
      guess.optns$type <- 'distance'
    }else{
      guess.optns$type <- 'compare'
    }

  }

  # fill-in if not specified
  if(is.null(optns$type)){
    optns$type <- guess.optns$type
  }
  optns$type <- match.arg(optns$type, c("distance", "thresholding", "compare"))

  # fill-in some default if not specified
  if(is.null(optns$family) & optns$type != "distance"){
    # if family not specified for thresholding and compare type
    optns$family <- stats::binomial(link = "logit")
  }

  if(is.null(optns$intercept)){
    optns$intercept <- base::ifelse(
      optns$type == "thresholding", TRUE, FALSE
    )
  }
  stopifnot(inherits(optns$intercept, 'logical'))

  return(optns)

}


#' Fit model for local metric
#'
#' @description
#' "fit" a model for local metric estimation, essentially prepare data.
#'
#' @param formula a formula for model
#' @param data a data frame
#' @param coord either a matrix of coordinates or missing
#' @param optns control options
#'
#' @return a \code{metricModel} object, essentially a list wrapping data and
#' options. Its elements are
#' \describe{
#'   \item{\code{graph}}{A list of graph edges for faster indexing.}
#'   \item{\code{coord}}{matrix of coordinates of points, one row is one point.}
#'   \item{\code{resp}}{
#'     data frame of response, in the format of \code{resp} as demanded by
#'     \code{\link{estiMetric}}.
#'   }
#'   \item{\code{optns}}{
#'     list of options, unspecified will be filled by default or best guesses.
#'   }
#' }
#' @export
#'
#' @details
#'
#' For squared distance type, \code{y} is numeric vector, modeled by
#' \deqn{E(y_{ij} | D_{ij}) = g^{-1}(\beta_0 + D_{ij}^2);}
#' for thresholding (binary censoring) type, \code{y} is binary vector, modeled
#' by
#' \deqn{P(y_{ij} = 1 | D_{ij}) = g^{-1}(c (\beta_0 + D_{ij}^2) );}
#' for comparing edges, \code{y} is binary vector, modeled by
#' \deqn{
#'     P(y_{ijkl} = 1 | D_{ij}, D_{kl})
#'     = g^{-1}(c (\beta_0 + D_{ij}^2 - D_{kl}^2) );
#' }
#' where
#' \eqn{y_{ij}} and \eqn{y_{ijkl}} are the response,
#' \eqn{D_{ij}} and \eqn{D_{kl}} are geodesic distance between i-j and k-l
#' points, while
#' \eqn{g} is the link function specified by \code{optns$family} (default to
#' Gaussian if missing),  \eqn{c} is some constant function due to conformality,
#' and \eqn{\beta_0} is intercept (optional). Note that here \eqn{c} cannot be
#' estimated, and the resulting estimation is conformal to the underlying truth.
#'
#' The model is specified by \code{formula}, use \cr
#' \code{y ~ (p1_1 + ... p1_d):(p2_1 + ... p2_d)} \cr
#' for squared distance or binary censoring responses;
#' use \cr\code{
#' y ~ (p1_1 + ... p1_d):(p2_1 + ... p2_d):(p3_1 + ... p3_d):(p4_1 + ... p4_d)
#' }\cr
#' for comparing edges with binary responses;
#' where the \code{y} are response variable in the \code{data}, while
#' \code{p1_1, ..., p4_d} are coordinates or edge indices in the
#' \code{data}, where the variables sandwiched between \code{:} defines a
#' point, \code{:} defines an edge or comparison of two edges. A group of two
#' consecutive points register one edge. For example, the previous formula
#' compares edge \code{p1-p2} with \code{p3-p4}. Note that the coordinates must
#' be in a same order for all the points.
#' \cr
#' As previously mentioned, one can use coordinates or indices in the bracket to
#' define a point. When using indices, coordinates of the points must be
#' supplied via \code{coord} as a matrix with one row being one point.
#'
#' The data frame \code{data} should contain response and edges, for example,
#' with formula \code{y ~ p1 : p2}, then in \code{data},
#' the \code{y} column is response, while \code{p1, p2} columns are index of
#' pairs of points, in which case, for example, a row of (100, 5, 6) means the
#' observed response between the 5th and 6th points (coordinates in 5th and 6th
#' row of \code{coord} is 100.
#' \cr
#' In addition to those described in \code{\link{set_optns}} and
#' \code{\link{guess_model}}, one additional control option is
#' \describe{
#'    \item{\code{local.reach}}{
#'        positive numbers defining range of local neighborhood near
#'        target points. If missing, will be set to approximately include 1000
#'        points in local neighbor hood.
#'    }
#' }
#'
#' @family {locMetric}
#'
#' @examples
#' d <- 3
#' manifold <- spaceEuclidean(d)
#' set.seed(1)
#' obsv.coord <- manifold$genPnt(10^4)
#' idx.edge <- allEdge(obsv.coord, local.reach = 0.1)
#' idx.edge <- idx.edge[idx.edge[, 1] != idx.edge[, 2], , drop = F]
#' idx.edge <- idx.edge[sample.int(nrow(idx.edge), 10^5), , drop = F]
#' resp <- manifold$dist(
#'   obsv.coord[idx.edge[, 1], ], obsv.coord[idx.edge[, 2], ]
#' ) ^ 2
#'
#' # input coord in data
#' data.w.coord <- cbind(
#'   resp, obsv.coord[idx.edge[, 1], ], obsv.coord[idx.edge[, 2], ]
#' )
#' data.w.coord <- as.data.frame(data.w.coord)
#' names(data.w.coord) <- c('y', sprintf('start%s', seq(d)), sprintf('end%s', seq(d)))
#' formula.w.coord <- as.formula(sprintf(
#'   "y ~ (%s) : (%s)",
#'   paste(sprintf('start%s', seq(d)), collapse = ' + '),
#'   paste(sprintf('end%s', seq(d)), collapse = ' + ')
#' ))
#'
#' # or use indices
#' data.w.idx <- cbind(resp, idx.edge)
#' data.w.idx <- as.data.frame(data.w.idx)
#' names(data.w.idx) <- c('y', sprintf('p%s', seq(2)))
#' formula.w.idx <- y ~ p1 : p2
#'
#' fit.w.coord <- fitMetric(formula.w.coord, data.w.coord)
#' fit.w.idx <- fitMetric(formula.w.idx, data.w.idx, coord = obsv.coord)
#'
#' all.equal(
#'   estiMetric(rep(0.5, d), fit.w.coord), estiMetric(rep(0.5, d), fit.w.idx)
#' )
#'
#' set.seed(10)
#' target <- matrix(runif(d * 10), ncol = d)
#' true.metric <- apply(target, 1, manifold$metric, simplify = F)
#'
#' # interchangeable model
#' all.equal(
#'   estiMetric(target, fit.w.coord), estiMetric(target, fit.w.idx)
#' )
#'
#' esti.metric <- estiMetric(target, fit.w.coord)
#' names(esti.metric) <- NULL
#'
#' # check symmetric and positive definite
#' all(sapply(esti.metric, isSymmetric))
#' all(sapply(esti.metric, function(metric){
#'   eig.val <- eigen(metric, symmetric = T, only.values = T)
#'   all(eig.val$values >= 0)
#' }))
#' # check estimation accuracy
#' all.equal(esti.metric, true.metric)
fitMetric <- function(formula, data, coord, optns = list()){
  # "fit" a model for local metric estimation, essentially prepare data.

  # # DEV ONLY
  # env.dev <- new.env()
  # local(envir = env.dev, {
  #
  #   # gen data
  #   d <- 3
  #   manifold <- spaceEuclidean(d)
  #   set.seed(1)
  #   obsv.coord <- manifold$genPnt(10^4)
  #   idx.edge <- allEdge(obsv.coord, local.reach = 0.1)
  #   idx.edge <- idx.edge[idx.edge[, 1] != idx.edge[, 2], , drop = F]
  #   idx.edge <- idx.edge[sample.int(nrow(idx.edge), 10^5), , drop = F]
  #   resp <- manifold$dist(
  #     obsv.coord[idx.edge[, 1], ], obsv.coord[idx.edge[, 2], ]
  #   ) ^ 2
  #   optns <- list()
  #   browser();QWER
  #   # input coord in data
  #   data <- cbind(
  #     resp, obsv.coord[idx.edge[, 1], ], obsv.coord[idx.edge[, 2], ]
  #   )
  #   data <- as.data.frame(data)
  #   names(data) <- c('y', sprintf('start%s', seq(d)), sprintf('end%s', seq(d)))
  #   formula <- as.formula(sprintf(
  #     "y ~ (%s) : (%s)",
  #     paste(sprintf('start%s', seq(d)), collapse = ' + '),
  #     paste(sprintf('end%s', seq(d)), collapse = ' + ')
  #   ))
  #
  #   # or use indices
  #   data <- cbind(resp, idx.edge)
  #   data <- as.data.frame(data)
  #   names(data) <- c('y', sprintf('p%s', seq(2)))
  #   formula <- y ~ p1 : p2
  #   coord <- obsv.coord
  # })
  # local(envir = env.dev, {browser();QWER})
  # # END DEV

  stopifnot(all(is.finite(as.matrix(data))))

  char.form <- deparse(formula)
  char.form <- stringr::str_remove_all(char.form, '\\s')
  # form.coord <- stringr::str_extract_all(
  #   char.form, "((?<=(:\\()).*(?=\\)))|((?<=(\\()).*(?=\\):))"
  # )
  form.coord <- stringr::str_remove_all(char.form, '.*~')
  form.coord <- stringr::str_split(form.coord, ':')
  form.coord <- sapply(unlist(form.coord), function(x) {
    stringr::str_split(stringr::str_remove_all(x, '\\(|\\)'), pattern = '\\+')
  })
  form.resp <- stringr::str_extract(char.form, '^.*(?=~)')

  # dimension
  if(!missing(coord)){
    d <- ncol(coord)
  }else{
    d <- length(form.coord[[1]])
    stopifnot(all(sapply(form.coord, length) == d))
  }

  # extract coordinates of input points
  input.pnts <- lapply(form.coord, function(i) {
    mat <- as.matrix(data[, i])
    colnames(mat) <- NULL
    return(mat)
  })

  # prepare mat.edge and obsv.coord
  if(!missing(coord)){

    mat.edge <- do.call(cbind, input.pnts)
    # some checking
    stopifnot(all(mat.edge >= 1) & all(mat.edge <= nrow(coord)))
    stopifnot(all.equal(as.integer(mat.edge), as.numeric(mat.edge)))

    stopifnot(all(is.finite(coord)))
    obsv.coord <- coord

  }else{

    # create look-up dictionary as a whole
    tm <- indexPoints(do.call(rbind, input.pnts))
    # this is the matrix of unique coordinates, ncol(coord) = d
    obsv.coord <- tm$coord
    idx.pnts <- tm$idx
    # split indices and reformat in to indices of edge
    idx.edge <- split(
      idx.pnts, rep(seq_along(input.pnts), each = nrow(data))
    )
    # now a matrix of edges, nrow(mat.edge) = nrow(data)
    mat.edge <- do.call(cbind, idx.edge)
    rm(idx.edge) # reserve name to avoid confusion.

  }
  resp <- cbind(data[[form.resp]], mat.edge)

  ### check options, keep identical as those in estiMetric
  if(is.null(optns$local.reach)){
    # approx. 1000 pnts in each nbhd, if target & obsv unif in a square
    optns$local.reach <-
      abs( prod( apply( obsv.coord, 2, function(x) { diff(range(x)) }) ) )  /
      (nrow(obsv.coord) / 1000)
  }
  optns <- set_optns(optns)
  optns <- guess_model(optns, formula = formula, data = resp)

  # if(is.null(optns$n.local)){
  #   optns$n.local <- c(10, 1000)
  # }
  # if(length(optns$n.local) == 1){
  #   optns$n.local <- c(1, optns$n.local)
  # }
  # stopifnot(all(optns$n.local >= 1))
  # if(is.null(optns$get.loc.obsv))
  #   optns$get.loc.obsv <- FALSE
  # stopifnot(
  #   length(optns$get.loc.obsv) == 1 & inherits(optns$get.loc.obsv, 'logical')
  # )
  # if(is.null(optns$locMetric)){
  #   tm <- list()
  #   if( length(unique(resp[, 1])) > 2 ){
  #     if(ncol(resp) != 3) stop(
  #       'input distance, but with edge comparison format.'
  #     )
  #     tm$type <- 'distance'
  #   }else if( ncol(resp) == 3 ){
  #     tm$type <- 'thresholding'
  #   }else{
  #     tm$type <- 'compare'
  #   }
  #   optns$locMetric <- tm
  # }
  # optns$method.trim <- match.arg(optns$method.trim, c("proximity", "random"))

  ### reformat for faster slicing later
  # for faster slicing later, we translate df.resp into list, for the record,
  # slicing list is 10 times faster than slicing data.frame, and slicing list
  # by its index is 100 times faster than slicing by its name. We compute diff
  # of coord here so as to avoid it being computed multiple times in different
  # local nbhd when estimating locally, especially local.reach is large creating
  # overlaps, and when comparing multiple pairs.

  df.resp <- data.frame(resp)
  names(df.resp) <- c('resp', outer(
    seq(2), seq(ncol(df.resp) %/% 2),
    function(x, y) sprintf('p%se%s', x, y)
  ))
  df.resp$idx.dist <- seq(nrow(df.resp))
  # next: to list, essentially a faster slicing index by split into list
  # this template of list is used for faster slicing later
  template.ls <- rep(list(NULL), nrow(obsv.coord))
  names(template.ls) <- seq_along(template.ls)

  ls.graph <- lapply(seq(ncol(resp) %/% 2), function(idx.edge.set) {
    # browser();QWER
    mat.edge <- as.matrix(df.resp[, 1 + seq(2) + 2 * (idx.edge.set - 1)])
    if(length(unique(mat.edge[, 1]))> length(unique(mat.edge[, 2]))){
      # switch index columns to put the least varying first
      mat.edge <- mat.edge[, c(2, 1)]
    }
    tm <- list(
      mat.edge, df.resp$idx.dist, df.resp$resp,
      deltaCoord(obsv.coord, mat.edge, return.format = 'keep')
    )
    tm <- do.call(cbind, tm) # this is fast enough
    colnames(tm) <- c(
      'p1', 'p2', 'idx.dist', 'resp', sprintf('coordd.%s', seq(d))
    )

    # split into list
    tm <- split.data.frame(tm, tm[, 'p1'])
    res <- template.ls
    res[names(tm)] <- tm
    return(res)
  })

  res <- list(
    graph = ls.graph, coord = obsv.coord,
    resp = df.resp, optns = optns
  )
  class(res) <- c('metricModel', class(res))

  return(res)

}

#' Local Estimation of Metric Tensor
#' @description
#'   Local estimation of metric tensor given squared distance,
#'   binary censoring, or comparing edges responses. As a high-level wrapper
#'   for \code{locMetric}.
#'
#' @param target matrix of the coordinate of target points (1 row = 1 point)
#' @param model a fitted model from \code{fitMetric}
#' @param resp data.frame or matrix of 3 or 5 columns with rows of observations
#' @param obsv.coord coordinate matrix of observed points, 1 row = 1 point
#' @param optns a list of control options.
#'
#' @return a list, \code{"metric"} containing estimated metric tensor matrices
#' and \code{"optns"}, formatted local observations (\code{"loc.obsv"}) will be
#' provided if requested.
#'
#' @details
#' Call by either \cr
#' \code{estiMetric(target, model, optns = list())}
#' \cr or \cr
#' \code{estiMetric(target, resp, obsv.coord, optns = list())}.
#' \cr\cr
#' When using the second way to call, use the following input format.
#' For \code{resp},
#'    its first column is response, while 2nd to last columns are index of pairs
#'    of points. For example, a row of (100, 5, 6) means the observed
#'    response between the 5th and 6th points (coordinates in 5th and 6th row of
#'    \code{obsv.coord} is 100; a row of (1, 5, 6, 6, 8) means (only in the
#'    comparing edges cases) the response for comparing edge 5-6 and edge 6-8 is
#'    1.
#'    \cr
#'    If requested by setting \code{optns$metric.only} as \code{FALSE}, the
#'    returned value have an additional elements for local observations (edges),
#'    which is also a named list named after row indices of \code{target}, whose
#'    elements are matrices with columns corresponding to indices of edges
#'    endpoints, response, and difference in coordinates.
#'    \cr\cr
#' In addition to those described in \code{\link{set_optns}} and
#' \code{\link{guess_model}}, one additional control option is
#' \describe{
#'    \item{\code{local.reach}}{
#'        positive numbers defining range of local neighborhood near
#'        target points. If missing, will be set to approximately include 1000
#'        points in local neighbor hood.
#'    }
#' }
#'
#' @export
#'
#' @family {locMetric}
#' @examples TBD
estiMetric <- function(target, model, resp, obsv.coord, optns = list()){
  # function estimating metric tensor at target points
  # args:
  #   target: coordinate matrix of the target points (1 row = 1 point)
  #   obsv.coord: coordinate matrix of observed points, 1 row = 1 point
  #   resp: a data.frame or matrix of 3 or 5 columns with 1 row for 1
  #     observed pair where 2nd--last columns are index of pairs of points and
  #     the first is response, for example, a row of (100, 5, 6) means the
  #     observed response between the 5th and 6th points (coordinates in 5th
  #     and 6th row of obsv.coord) is 100; a row of (1, 5, 6, 6, 8) means (only
  #     in the comparing edges cases) the response for comparing edge 5-6 and
  #     edge 6-8 is 1.
  #     Note that the 2nd column should always be no smaller than the 1st.
  #     (TBD: Duplication will be checked. )
  #   optns: a list of control options.
  #     local.reach: maximum difference in coordinates for obsv.coord to target
  #       to still be consider in local neighbourhood of the target point.
  #     n.local: min and max number of points to be included in local nbhd.
  #     method.trim: "random"/"proximity", how to select local edges near target
  #     get.loc.obsv: whether to return local observations used
  #     locMetric: list passed to optns argument in locMetric function
  #     (TBD: check11: FALSE/TRUE(default), check resp for duplication.)
  # return: a list of estimated metric tensor matrices

  # DEBUG ONLY
  # browser();
  # env.dev <- new.env()
  # local(envir = env.dev, {
  #   optns <- list()
  #   d <- 3
  #   target <- rep(0, d)
  #   manifold <- spaceEuclidean(d)
  #   set.seed(1)
  #   obsv.coord <- manifold$genPnt(10000)
  #   idx.edge <- allEdge(obsv.coord, local.reach = 0.1)
  #   idx.edge <- idx.edge[idx.edge[, 1] != idx.edge[, 2], , drop = F]
  #   resp <- cbind(manifold$dist(
  #     obsv.coord[idx.edge[, 1], ], obsv.coord[idx.edge[, 2], ]
  #   ) ^ 2, idx.edge)
  #   browser();QWER
  #   plotGraph(obsv.coord, resp[, -1], max.n.edge = 50)
  # })
  # local(envir = env.dev, browser())
  # END DEBUG ONLY

  if(!missing(model)){
    if(!missing(resp) | !missing(obsv.coord))
      stop('Too many input, provide either 1. model or 2. resp and obsv.coord.')
  }else{

    # check coordinates
    stopifnot(is.matrix(obsv.coord))
    stopifnot(all(is.finite(obsv.coord)))

    # check distance
    stopifnot(ncol(resp) %in% c(3, 5))
    stopifnot(all(is.finite(resp)) & all(resp >= 0))

    # check indices
    stopifnot(all(
      resp[, -1] >= 1 & resp[, -1] <= nrow(obsv.coord)
    ))
    stopifnot(all.equal(
      as.numeric(resp[, -1]), as.integer(resp[, -1])
    ))

    # call fitMetric
    tm.dat <- data.frame(resp)
    names(tm.dat) <- c('y', sprintf('p%s', seq(ncol(resp) - 1)))
    model <- fitMetric(
      formula = as.formula(sprintf('y ~ %s', paste(
        sprintf('p%s', seq(ncol(resp) - 1)), collapse = ' : '
      ))),
      data = tm.dat, coord = obsv.coord, optns = optns
    )
    rm(tm.dat)

  }

  d <- ncol(model$coord)

  # check target
  if(!is.matrix(target)){
    stopifnot(length(target) == d)
    target <- matrix(target, nrow = 1)
  }else{
    stopifnot(ncol(target) == d)
  }

  ### check options
  optns <- model$optns
  if(is.null(optns$local.reach)){
    # approx. 1000 pnts in each nbhd, if target & obsv unif in a square
    optns$local.reach <-
      abs( prod( apply( model$coord, 2, function(x) { diff(range(x)) }) ) ) /
      (nrow(model$coord) / 1000)
  }
  optns <- set_optns(optns)
  optns <- guess_model(optns, formula = formula, data = resp)
  # if(is.null(optns$n.local)){
  #   optns$n.local <- c(10, 1000)
  # }
  # if(length(optns$n.local) == 1){
  #   optns$n.local <- c(1, optns$n.local)
  # }
  # stopifnot(all(optns$n.local >= 1))
  # if(is.null(optns$metric.only))
  #   optns$metric.only <- TRUE
  # stopifnot(
  #   length(optns$metric.only) == 1 & inherits(optns$metric.only, 'logical')
  # )
  # if(is.null(optns$locMetric)){
  #   tm <- list()
  #   if( length(unique(resp[, 1])) > 2 ){
  #     if(ncol(resp) != 3) stop(
  #       'input distance, but with edge comparison format.'
  #     )
  #     tm$type <- 'distance'
  #   }else if( ncol(resp) == 3 ){
  #     tm$type <- 'thresholding'
  #   }else{
  #     tm$type <- 'compare'
  #   }
  #   optns$locMetric <- tm
  # }
  optns$method.trim <- match.arg(optns$method.trim, c("proximity", "random"))

  local.reach <- optns$local.reach
  max.n.local <- max(optns$n.local)
  min.n.local <- min(optns$n.local)

  # df.resp <- data.frame(resp)
  # names(df.resp) <- c('resp', outer(
  #   seq(2), seq(ncol(df.resp) %/% 2),
  #   function(x, y) sprintf('p%se%s', x, y)
  # ))
  # df.resp$idx.dist <- seq(nrow(df.resp))
  # # next: to list, essentially a faster slicing index by split into list
  # # this template of list is used for faster slicing later
  # template.ls <- rep(list(NULL), nrow(obsv.coord))
  # names(template.ls) <- seq_along(template.ls)
  #
  # ls.graph <- lapply(seq(ncol(resp) %/% 2), function(idx.edge.set) {
  #   # browser();QWER
  #   mat.edge <- as.matrix(df.resp[, 1 + seq(2) + 2 * (idx.edge.set - 1)])
  #   if(length(unique(mat.edge[, 1]))> length(unique(mat.edge[, 2]))){
  #     # switch index columns to put the least varying first
  #     mat.edge <- mat.edge[, c(2, 1)]
  #   }
  #   tm <- list(
  #     mat.edge, df.resp$idx.dist, df.resp$resp,
  #     deltaCoord(obsv.coord, mat.edge, return.format = 'keep')
  #   )
  #   tm <- do.call(cbind, tm) # this is fast enough
  #   colnames(tm) <- c(
  #     'p1', 'p2', 'idx.dist', 'resp', sprintf('coordd.%s', seq(d))
  #   )
  #
  #   # split into list
  #   tm <- split.data.frame(tm, tm[, 'p1'])
  #   res <- template.ls
  #   res[names(tm)] <- tm
  #   return(res)
  # })
  # # object.size(ls.graph) / 1024 ^ 2

  ls.graph <- model$graph
  obsv.coord <- model$coord

  ### use diff in coord to determine local nbhd of target points
  target.nbhd <- locWindow(target, obsv.coord, local.reach)

  ## then two ways to subset edges: random or "proximity"
  if(optns$method.trim == 'proximity'){
    # pseudo distance to target points
    proximity <- lapply(seq_along(target.nbhd), function(idx.target){
      # browser();QWER
      tm <- obsv.coord[target.nbhd[[idx.target]], ]
      tm <- tm - target[idx.target, col(tm)]
      return(rowSums(tm ^ 2))
    })
    names(proximity) <- names(target.nbhd)
    target.nbhd <- base::mapply(function(nbhd, prox){
      # browser();QWER
      res <- nbhd[head(order(prox), max.n.local)] # now ordered by proximity
      res <- sort(res) # now by smaller indices to larger ones
      return(res)
    }, nbhd = target.nbhd, prox = proximity, SIMPLIFY = FALSE)
  }else{
    target.nbhd <- lapply(target.nbhd, function(x) {
      res <- sample(x, size = min(length(x), max.n.local), replace = F)
      return(sort(res))
    })
  }
  # slice the observations used for local estimation
  loc.obsv <- lapply(target.nbhd, function(nbhd){
    # browser();QWER

    loc.graph <- lapply(ls.graph, function(x) x[as.numeric(nbhd)])
    loc.graph <- lapply(loc.graph, function(x) {
      # browser();QWER
      res <- x[!sapply(x, is.null)]
      res <- lapply(res, function(mat) mat[mat[, 'p2'] %in% nbhd, , drop = F])
      res <- do.call(rbind, res)
      return(res)
    })
    if(length(loc.graph) == 2){
      # if comparing edges, there will be 2 lists, need to intersect
      idx.yes <- base::intersect(
        loc.graph[[1]][, 'idx.dist'], loc.graph[[2]][, 'idx.dist']
      )
      loc.graph <- lapply(loc.graph, function(mat){
        res <- mat[mat[, 'idx.dist'] %in% idx.yes, , drop = F]
        res <- res[order(res[, 'idx.dist']), , drop = F]
        res
      })
      return(loc.graph)
    }else{
      return(loc.graph[[1]])
    }

  })

  # check if densely sampled near targets
  if(length(ls.graph) == 2)
    n.loc.obsv <- sapply(loc.obsv, function(x) nrow(x[[1]]))
  else
    n.loc.obsv <- sapply(loc.obsv, nrow)
  if(!all(n.loc.obsv >= min.n.local)) {
    warning(sprintf(
      'samples are not dense enough near %s target points.',
      sum(n.loc.obsv < min.n.local)
    ))
  }

  # estimating for each targeted points
  res.locMetric <- lapply(loc.obsv, function(input){
    if(is.list(input)){
      # if comparing
      ls.args <- list(
        y = input[[1]][, 'resp'],
        coord.diff = lapply(input, function(x) x[, -seq(4), drop = F])
      )
    }else{
      ls.args <- list(
        y = input[, 'resp'],
        coord.diff = input[, -seq(4), drop = F]
      )
    }
    ls.args$optns <- optns
    return(do.call(locMetric, ls.args))
  })

  # return list of local results: estimated metric tensor only or not
  if(optns$metric.only)
    return(lapply(res.locMetric, `[[`, 'mat.g'))

  loc.res <- list(
    metric = lapply(res.locMetric, `[[`, 'mat.g'),
    obsv <- loc.obsv, optns = optns
  )
  return(loc.res)

}


#' Local Estimation of Metric Tensor
#' @description
#'   Local estimation of Riemannian metric tensor with noisy geodesic distances.
#' @param y
#'   array of response
#' @param coord.diff
#'   matrix of difference of coordinates corresponding to edges, one row is one
#'   edge (connecting a pair of points). Number of rows = length of \code{y}.
#'   If comparing edges, then a list of 2 such matrices.
#' @param optns
#'   a list of control options.
#'
#' @return a list of
#' \describe{
#'   \item{\code{mat.g}}{
#'     estimated metric tensor as ncol(coord.diff) X ncol(coord.diff) matrix.
#'    }
#'   \item{\code{mat.design}}{
#'     design matrix used for local likelihood.
#'    }
#'   \item{\code{loc.fit}}{
#'     glm.fit result.
#'    }
#'   \item{\code{optns}}{
#'     Options used.
#'    }
#' }
#' @details Possible control options are
#' \describe{
#'   \item{\code{type}}{
#'     \code{"distance"}(default), \code{"thresholding"}, or \code{"compare"},
#'     simplified input handled by \code{base::match.arg}.
#'     If \code{"compare"}, then \code{coord.diff} needs to be a list of two
#'     matrices each being difference of coordinates for the pair of points
#'     compared.
#'    }
#'   \item{\code{intercept}}{
#'     \code{FALSE} or \code{TRUE}, whether to add an intercept term.
#'     Will default to \code{FALSE} unless \code{type = "thresholding"}.
#'    }
#'   \item{\code{...}}{
#'     passed to \code{glm.fit}, if no \code{family} specified, will use a
#'     simple QR-solve for OLS. One can also specify weight when using
#'     \code{glm.fit}. Default:
#'     missing for \code{type = "distance"};
#'     \code{binomial(link = "logit")} for
#'        \code{type = c("thresholding", "compare")}.
#'    }
#' }
#' One can use \code{\link{guess_model}} to set those options.
#'
#' @export
#'
#' @family {locMetric}
#'
#' @examples
#' ## Euclidean space
#' set.seed(1)
#' # coord of points
#' coord.pnts <- matrix(runif(500 * d), ncol = d)
#' # pairing points
#' idx.dist <- t(utils::combn(nrow(coord.pnts), 2))
#' idx.dist <- idx.dist[
#'   sample(nrow(idx.dist), size = 5000),
#' ]
#' coord.diff <-
#'   coord.pnts[idx.dist[, 1], ] - coord.pnts[idx.dist[, 2], ]
#' arr.dist <- sqrt(rowSums(coord.diff ^ 2))
#' prob.thre <- sigmoid.f(arr.dist ^ 2)
#' prob.thre.itcpt <- sigmoid.f(arr.dist ^ 2 - mean(arr.dist ^ 2))
#' arr.thre <- rbinom(length(arr.dist), 1, prob = prob.thre)
#' arr.thre.itcpt <- rbinom(length(arr.dist), 1, prob = prob.thre.itcpt)
#'
#' # noiseless observation
#' locMetric(arr.dist ^ 2, coord.diff)$mat.g
#'
#' # normal noise (after squared)
#' locMetric(
#'   arr.dist ^ 2 + rnorm(length(arr.dist), sd = sd(arr.dist) / 50),
#'   coord.diff
#' )$mat.g
#'
#' # binary thresholding
#'  locMetric(
#'    arr.thre, coord.diff
#'    , optns = list(type = 't', intercept = FALSE)
#'  )$mat.g
#'
#' # comparing distance
#' set.seed(1)
#' idx.comp <- t(utils::combn(sample(length(arr.dist), 1000), 2))
#' idx.comp <- idx.comp[sample(nrow(idx.comp), size = 10000), ]
#' ls.coord.diff <- list(
#'   coord.diff[idx.comp[, 1], ],
#'   coord.diff[idx.comp[, 2], ]
#' )
#' prob.comp <- sigmoid.f(
#'   arr.dist[idx.comp[, 1]]^2 - arr.dist[idx.comp[, 2]]^2
#' )
#' arr.comp <- rbinom(length(prob.comp), 1, prob = prob.comp)
#'
#' # comparing distance
#' locMetric(
#'   arr.comp, ls.coord.diff
#'   , optns = list(type = 'comp', intercept = FALSE)
#' )$mat.g
#' # comparing distance w/ intercept
#' res.itcpt <- locMetric(
#'   arr.comp, ls.coord.diff
#'   , optns = list(type = 'comp', intercept = TRUE)
#' )
#' res.itcpt$mat.g
#' # the fitted intercept should be close to 0
#' res.itcpt$loc.fit$coefficients[1]
locMetric <- function(y, coord.diff, optns = list()){
  # local polynomial estimating metric tensor.
  # browser()

  # optns
  if(is.null(optns$type)){
    optns$type <- "distance"
  }
  optns$type <- match.arg(optns$type, c("distance", "thresholding", "compare"))
  if(is.null(optns$family) & optns$type != "distance"){
    # if family not specified for thresholding and compare type
    optns$family <- stats::binomial(link = "logit")
  }

  if(is.null(optns$intercept)){
    optns$intercept <- base::ifelse(
      optns$type == "thresholding", TRUE, FALSE
    )
  }
  # idx of arguments not pass to glm.fit
  # idx.optns.notglm <- which(names(optns) %in% c('intercept', 'type'))
  idx.optns.glm <- idx_glm_optns(optns)

  # sanity and desgin matrix
  if(optns$type != "compare"){
    stopifnot(is.matrix(coord.diff))
    stopifnot(length(y) == nrow(coord.diff))

    coord.d <- ncol(coord.diff)
    mat.v <- vecQuad(coord.diff)
  }else{
    stopifnot(is.list(coord.diff))
    stopifnot(length(coord.diff) == 2)
    stopifnot(is.matrix(coord.diff[[1]]) & is.matrix(coord.diff[[2]]))
    stopifnot(all(dim(coord.diff[[1]]) == dim(coord.diff[[2]])))
    stopifnot(length(y) == nrow(coord.diff[[1]]))

    coord.d <- ncol(coord.diff[[1]])
    mat.v <- vecQuad(coord.diff[[1]]) - vecQuad(coord.diff[[2]])
  }

  if(optns$intercept == TRUE){
    mat.v <- cbind(rep(1, nrow(mat.v)), mat.v) # add intercept if asked to
  }

  if(!is.null(optns$family)){# if using glm.fit with general family
    loc.fit <- do.call(stats::glm.fit, c(
      list(x = mat.v, y = y),
      optns[idx.optns.glm]
    ))
    vec.g <- loc.fit$coefficients
  }else{
    # vec.g <- solve(crossprod(mat.v)) %*% colSums(mat.v * y)
    vec.g <- qr.solve(crossprod(mat.v), colSums(mat.v * y))
    loc.fit <- NULL
  }
  if(optns$intercept == TRUE){
    vec.g <- vec.g[-1] # remove fitted intercept from estimated metric
  }
  # translate to metric tensor
  mat.g <- matrix(0, coord.d, coord.d)
  mat.g[lower.tri(mat.g, diag = FALSE)] <- #lower.tri since symmetric
    vec.g[-seq(coord.d)] / 2
  mat.g <- t(mat.g) + mat.g
  diag(mat.g) <- vec.g[seq(coord.d)]

  return(list(
    mat.g = mat.g,
    mat.design = mat.v,
    loc.fit = loc.fit
    , optns = optns
  ))
}
