# core functions for metric estimation

#' Local Estimation of Metric Tensor
#' @description
#'   Local estimation of Riemannian metric tensor with noisy geodesic distances.
#' @param y
#'   array of response
#' @param coord.diff
#'   difference of coordinates, one row is one pair of points.
#'   Number of rows equal to length of \code{y}.
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
#'     passed to \code{glm.fit}, if non, will use a simple QR-solve for OLS.
#'     Default: missing for \code{type = "distance"};
#'       \code{binomial(link = "logit")} for
#'       \code{type = c("thresholding", "compare")}.
#'    }
#' }
#'
#' For distance type, \code{y} is numeric vector, modeled by
#' \deqn{E(y_i | D_i) = g^{-1}(\beta_0 + D_i^2),}
#' for thresholding type, \code{y} is binary vector, modeled by
#' \deqn{P(y_i = 1 | D_i) = g^{-1}(\beta_0 + D_i^2),}
#' for compare type, \code{y} is binary vector, modeled by
#' \deqn{P(y_{ij} = 1 | D_i) = g^{-1}(\beta_0 + D_i^2 - D_j^2),}
#' where
#' \eqn{y_i} is the ith response, \eqn{D_i} is the ith pair of geodesic distance
#' corresponding to ith row of \code{coord.diff},
#' \eqn{g} is the link function specified by \code{optns$family} (default to
#' Gaussian if missing), \eqn{\beta_0} is intercept (optional).
#' @export
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
#' locMetric.g(arr.dist ^ 2, coord.diff)$mat.g
#'
#' # normal noise (after squared)
#' locMetric.g(
#'   arr.dist ^ 2 + rnorm(length(arr.dist), sd = sd(arr.dist) / 50),
#'   coord.diff
#' )$mat.g
#'
#' # binary thresholding
#'  locMetric.g(
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
#' locMetric.g(
#'   arr.comp, ls.coord.diff
#'   , optns = list(type = 'comp', intercept = FALSE)
#' )$mat.g
#' # comparing distance w/ intercept
#' res.itcpt <- locMetric.g(
#'   arr.comp, ls.coord.diff
#'   , optns = list(type = 'comp', intercept = TRUE)
#' )
#' res.itcpt$mat.g
#' # the fitted intercept should be close to 0
#' res.itcpt$loc.fit$coefficients[1]
locMetric.g <- function(y, coord.diff, optns = list()){
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
  idx.optns.notglm <- which(names(optns) %in% c('intercept', 'type'))

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
    loc.fit <- do.call(glm.fit, c(
      list(x = mat.v, y = y),
      optns[-idx.optns.notglm]
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
