### functions for local likelihood

### parts: x is the target point
###   1. get model matrix,
###     input: obsv, local poly specification, output: a func. of (x)
###   2. the weights:
###     a. kernel specification: get K_h() from kernel(K) and bandwidth(h);
###     b. get weight matrix, take output of a. and obsv, output a func.(x)
###   3. objective function constructor: (only needed for quasi likelihood?)
###     input: "link" function: loc.poly <--> likelihood
###     output: a func.(x, beta), where beta are parameters to optimize.

#' Define a local polynomial.
#'
#' @param formula a formula for local polynomial.
#' @param data a data frame (or matrix, if terms are column number).
#' @param optns a list of control options.
#'
#' @details
#' The local polynomial will be defined by the right-hand-side of the
#' \code{formula}, for example, the corresponding model matrix for
#' \code{y ~ a * b} at every local point \code{x} is
#' \code{model.matrix(~(a-x) * (b-y))} following default interpretation of R.
#'
#' Possible control options are
#' \describe{
#' \item{\code{kernel.f}}{a symmetric nonnegative univariate function, e.g.,
#' output of getKernelFunc.}
#'
#' \item{\code{zero.weight.eps}}{a small positive value s.t. weights smaller
#'   than it will be designated as zero and corresponding observation will be
#'   exclude in local estimation.}
#' }
#'
#' @return a list containing
#' \describe{
#'   \item{\code{mdl.mat.f}}{A function taking argument \code{(x, wt)} computes
#'   local model matrix at \code{x} with weight \code{wt}(optional).}
#'   \item{\code{weight.f}}{
#'     (Optional) If \code{optns$kernel.f} is provided, the function determine
#'     weight of observations at any given local point.
#'   }
#'   \item{\code{mdl.comp}}{A list of model components, \code{formula} is the
#'   defining formula, while \code{df.mdl} contains the part of input data
#'   involved. The \code{optns} is also recorded.}
#' }
#' @examples TBD
#' @export
locPolyDefine <- function(formula, data, optns = list()){
  # # DEBUG
  # formula <- y ~ a * b
  # data <- data.frame(
  #   y = seq(10), a = seq(10) + 100, b = seq(10) - 100
  #   , d = seq(10) * 50
  # )
  # browser();QWER
  # # END DEBUG

  stopifnot(is.matrix(data) | is.data.frame(data))
  data <- as.data.frame(data)

  # select the data according to terms
  o.mat <- as.matrix(stats::get_all_vars(formula[-2], data))
  df.resp <- stats::get_all_vars(formula[-3], data)
  # stopifnot(!anyNA(o.mat)) # no NA allowed.

  n.terms <- ncol(o.mat)
  n.obsv <- nrow(o.mat)

  mdl.mat.f <- function(x, wt){
    if(!is.null(names(x))){ # match by names if provided
      x <- x[colnames(o.mat)]
    }
    stopifnot(length(x) == n.terms)
    x <- as.numeric(x)
    w.df <- data.frame(df.resp[[1]], t(t(o.mat) - x))
    names(w.df) <- c(names(df.resp), colnames(o.mat))
    if(!missing(wt)){ # drop zero weight observations
      stopifnot(all(wt >= 0) & length(wt) == n.obsv)
      w.df <- w.df[wt != 0, , drop = FALSE]
    }
    res.mat <- stats::model.matrix(formula, w.df)
    return(res.mat)
  }

  # if specifying weight here
  if(is.null(optns$kernel.f)){
    weight.f <- NULL
  }else{
    if(is.null(optns$zero.weight.eps)){
      zero.weight.eps <- 0
    }else{
      zero.weight.eps <- optns$zero.weight.eps
    }
    stopifnot(length(zero.weight.eps) == 1 & zero.weight.eps >= 0)
    optns$zero.weight.eps <- zero.weight.eps
    weight.f <- function(x){
      # browser()
      if(!is.null(names(x))){ # match by names if provided
        x <- x[colnames(o.mat)]
      }
      stopifnot(length(x) == n.terms)
      x <- as.numeric(x)
      w.mat <- t(t(o.mat) - x)
      # L2 distance (TBD: anisotropy bandwidth)
      w.arr <- sqrt(rowSums(w.mat ^ 2))
      w.arr <- optns$kernel.f(w.arr)
      w.arr[w.arr < zero.weight.eps] <- 0
      return(w.arr)
    }
  }
  return(list(
    mdl.mat.f = mdl.mat.f,
    weight.f = weight.f,
    mdl.comp = list(
      df.mdl = cbind(df.resp, o.mat), formula = formula
      , optns = optns
    )
  ))
}


#' Get kernel function from name.
#'
#' @param kernel name of the kernel function, "gauss" or "step".
#' @param bandwidth bandwidth to use.
#' @param optns a list of control options (TBD).
#'
#' @return a univariate function \eqn{x \mapsto K(x/h)/h} where
#'   \eqn{h} is the bandwidth.
#'
#' @details For Gaussian kernel, the \eqn{K} is the Gaussian density function,
#'   for step kernel, the \eqn{K} is the step function on \eqn{[-1, 1]} taking
#'   value of 0.5.
#'
#' @examples
#' plot(getKernelFunc('gauss', bandwidth = 1, n.terms = 1), xlim = c(-5, 5))
#' plot(dnorm, xlim = c(-5, 5), col = 'red', add = TRUE)
#' plot(getKernelFunc('step', bandwidth = 5, n.terms = 1), xlim = c(-5, 5))
#' plot(dnorm, sd = 5, xlim = c(-5, 5), col = 'red', add = TRUE)
#' @export
getKernelFunc <- function(kernel = 'gauss', bandwidth = 1, optns = list()){
  stopifnot(bandwidth > 0)
  kernel <- match.arg(kernel, c("gaussian", 'step'))
  if(kernel == 'gaussian'){
    o.kf <- dnorm
  }
  if(kernel == 'step'){
    o.kf <- function(z) 0.5 * as.numeric(abs(z) <= 1)
  }
  return(function(x) o.kf(x / bandwidth) / bandwidth)
}
# plot(getKernelFunc('gauss', bandwidth = 1, n.terms = 1), xlim = c(-5, 5))
# plot(dnorm, xlim = c(-5, 5), col = 'red', add = TRUE)
# plot(getKernelFunc('step', bandwidth = 5, n.terms = 1), xlim = c(-5, 5))
# plot(dnorm, sd = 5, xlim = c(-5, 5), col = 'red', add = TRUE)

#' Get function computing weights
#'
#' @param loc.poly output of \code{\link[remeloc]{locPolyDefine}}.
#' @param kernel.f a symmetric nonnegative univariate function, e.g., output
#     of getKernelFunc.
#' @param optns a list of control options.
#'
#' @return a function computing weight as \code{kernel.f(obsv - x)}.
#'
#' @details Possible options are
#' \describe{
#' \item{\code{zero.weight.eps}}{a small positive value s.t. weights smaller
#'   than it will be designated as zero and corresponding observation will be
#'   exclude in local estimation.}
#' }
#' @examples TBD
#' @export
locPolyWeight <- function(loc.poly, kernel.f, optns = list()){

  n.terms <- ncol(loc.poly$mdl.comp$df.mdl) - 1
  o.mat <- as.matrix(loc.poly$mdl.comp$df.mdl[, -1])

  if(is.null(optns$zero.weight.eps)){
    zero.weight.eps <- 0
  }else{
    zero.weight.eps <- optns$zero.weight.eps
  }
  stopifnot(length(zero.weight.eps) == 1 & zero.weight.eps >= 0)

  res.f <- function(x){
    # browser()
    x <- as.numeric(x)
    stopifnot(length(x) == n.terms)
    w.mat <- t(t(o.mat) - x)
    # L2 distance (TBD: anisotropy bandwidth)
    w.arr <- sqrt(rowSums(w.mat ^ 2))
    w.arr <- kernel.f(w.arr)
    w.arr[w.arr < zero.weight.eps] <- 0
    return(w.arr)
  }
  return(res.f)
}

### DEVELOPING

# TBD, maybe not necessary.
# locPolyEsti <- function(loc.poly, newdata, optns = list()){
#   # local polynomial estimation
#   # args:
#   #   loc.poly: output of locPolyDefine
#   #   newdata: newdata
#   #   optns: list of options
#   stopifnot(all(loc.poly$terms %in% colnames(newdata)))
#
#   w.df <- as.data.frame(newdata[loc.poly$terms])
#
# }

# tst <- factor(x = c(3,2,4,5), levels = seq(6))
# tst <- factor(x = c('a', 'd'), levels = c('d', 'b', 'a'))
# tst
# tst %>% as.integer
# levels(tst)
multiNomialGLM <- function(y, x, wt){
  # fitting GLM for catagorical/multinomial family with logistic link
  # args:
  #   y: a factor/character vector, 1 element is one response.
  #   x: covariate matrix, same nrow as length(y). If want intercept, put col of
  #      1s in here.
  #   wt: weight for each reponse.
  # returns:
  #   fitted coefficients.
  # Details:
  #   P(Y = j | X) = \pi_j with \pi_j = exp(\eta_j) / \sum_l exp(\eta_l)
  #   where \eta_j = linear combination of X.

  # # DBUG
  # y <- factor(
  #   # x = sample(c('a', 'b', 'c'), size = 10, replace = TRUE)
  #   x = c('a', 'a', 'a', 'a', 'a', 'b', 'b', 'c', 'c', 'c')
  #   , levels = c('b', 'a', 'c')
  #   # , levels = c('c', 'a')
  # )
  # set.seed(1)
  # # x <- matrix(c(rep(1, length(y)), runif(length(y))), ncol = 2)
  # x <- matrix(c(rep(1, length(y)), seq(length(y))), ncol = 2)
  # # eta <-
  # # sample(
  # #   c('a', 'b'), size = nrow(x), replace = TRUE
  # #   , prob =
  # # )
  # # END DEBUG

  if(missing(wt)) wt <- rep(1, length(y))
  stopifnot(!(anyNA(y) | anyNA(x) | anyNA(wt)))
  stopifnot(is.matrix(x))
  stopifnot((length(y) == nrow(x)) & (length(y) == length(wt)))
  y <- as.factor(y)

  # unused levels will be dropped in computation.
  input.y <- y
  unseen.lvl <- setdiff(levels(y), unique(y))
  if(length(unseen.lvl) != 0){
    warning(
      sprintf(
        'Levels %s unseen, fitted coefficients could be unreliable.'
        , paste0(unseen.lvl, collapse = ', ')
      )
    )
    y <- base::droplevels(y)
  }
  n.lvl <- nlevels(y)
  stopifnot(n.lvl >= 1)
  idx.input.y.lvl <- base::match(levels(y), table = levels(input.y))
  # index for response in levels
  idx.y.lvl <- as.integer(y)

  if(n.lvl == 1){ # only observe 1 category
    warning('Only observing 1 level in y.')

    mat.beta <- matrix(0, nrow = ncol(x), ncol = 1)
    mat.eta <- matrix(0, nrow = nrow(x), ncol = 1)
    mat.pi <- matrix(1, nrow = nrow(x), ncol = 1)
    optim.res <- NULL
  }else{ # more than 1 categories

    # matrix for coefficients, one column for 1 category
    mat.beta <- matrix(0, nrow = ncol(x), ncol = n.lvl - 1)
    # # matrix for predictor eta, we set eta_J = 0 for identifiability
    # mat.eta <- matrix(0, nrow = nrow(x), ncol = n.lvl)
    # mat.eta[, seq(n.lvl - 1)] <- tcrossprod(x, t(mat.beta))
    # # matrix for prob. of each category, 1 row for 1 response, not really needed
    # mat.pi <- exp(mat.eta)
    # mat.pi <- mat.pi / rowSums(mat.pi)

    if(ncol(x) == 1 & all(base::duplicated(x)[-1])){ # intercept only
      # linear predictor is just weighted average, take log to get coefficients
      # though recall we set the eta (now = intercept) for last category as 0.
      beta.last <- log(sum(wt[idx.y.lvl == n.lvl]))
      for(idx.lvl in seq(n.lvl - 1)){
        mat.beta[, idx.lvl] <- log(sum(wt[idx.y.lvl == idx.lvl])) - beta.last
      }
      optim.res <- NULL
    }else{
      # for indexing mat.eta s.t. t(mat.eta)[arr.idx.eta[i] + j] = mat.eta[i, j]
      arr.idx.eta <- seq(0, n.lvl * (nrow(x) - 1), by = n.lvl)

      # negativ log-likelihood function
      neglogll <- function(vec.beta){
        mat.beta <- matrix(vec.beta, nrow = ncol(x), ncol = n.lvl - 1)
        mat.eta <- matrix(0, nrow = nrow(x), ncol = n.lvl)
        mat.eta[, seq(n.lvl - 1)] <- tcrossprod(x, t(mat.beta))
        arr.eta <- t(mat.eta)[arr.idx.eta + idx.y.lvl]
        return(
          -1 * sum(wt * (arr.eta - log(rowSums(exp(mat.eta)))))
        )
      }
      optim.res <- optim(as.numeric(mat.beta), neglogll)#, method = 'L-BFGS-B')

      # fitted coefficients
      mat.beta[, seq(n.lvl - 1)] <- optim.res$par
    }

    # corresponding predictor eta and probabilities
    mat.eta <- matrix(0, nrow = nrow(x), ncol = n.lvl)
    mat.eta[, seq(n.lvl - 1)] <- tcrossprod(x, t(mat.beta))
    mat.pi <- exp(mat.eta)
    mat.pi <- mat.pi / rowSums(mat.pi)
    mat.beta <- cbind(mat.beta, 0) # the zeros previously dropped for identifiable
  }

  # naming
  colnames(mat.beta) <- levels(y)
  colnames(mat.pi) <- levels(y)
  colnames(mat.eta) <- levels(y)
  rownames(mat.beta) <- colnames(x)

  # put back the unused levels
  mat.beta.full <- matrix(-Inf, nrow = ncol(x), ncol = nlevels(input.y))
  mat.beta.full[, idx.input.y.lvl] <- mat.beta
  mat.eta.full <- matrix(-Inf, nrow = nrow(x), ncol = nlevels(input.y))
  mat.eta.full[, idx.input.y.lvl] <- mat.eta
  mat.pi.full <- matrix(0, nrow = nrow(x), ncol = nlevels(input.y))
  mat.pi.full[, idx.input.y.lvl] <- mat.pi
  # naming
  colnames(mat.beta.full) <- levels(input.y)
  colnames(mat.pi.full) <- levels(input.y)
  colnames(mat.eta.full) <- levels(input.y)
  rownames(mat.beta.full) <- colnames(x)

  # y
  # levels(y)[mat.pi %>% apply(1, which.max)]
  # factor(levels(y)[apply(mat.pi, 1, which.max)], levels = levels(y))

  return(list(
    coefficients = mat.beta.full,
    fitted.values = mat.pi.full,
    fitted.y =
      factor(levels(y)[apply(mat.pi, 1, which.max)], levels = levels(y)),
    linear.predictors = mat.eta.full,
    optim.res = optim.res
  ))
}
# y <- factor(
#   # x = sample(c('a', 'b', 'c'), size = 10, replace = TRUE)
#   x = c('a', 'a', 'a', 'a', 'a', 'b', 'b', 'c', 'c', 'c')
#   , levels = c('b', 'a', 'c', 'd')
#   # , levels = c('c', 'a')
# )
# x <- matrix(c(rep(1, length(y)), seq(length(y))), ncol = 2)
# multiNomialGLM(y, x)

locPoly_getModelInfo <- function(family = 'multinomial'){
  # function wrapper for model list to use with caret package.

  family <- match.arg(family, c(
    'multinomial' # TBD, 'gaussian'
  ))

  res <- list(
    library = NULL, type = 'Classification'
    , parameters = data.frame(
      parameter = c('poly.deg', 'kernel.bandwidth'),
      class = c('numeric', 'numeric'),
      label = c('polynomial degree', 'kernel bandwidth')
    )
    , grid = function(x, y, len = NULL, search = 'grid'){
      # a function that is used to create the tuning grid
      # n.terms <- ncol(x)
      mid.diff <- mean(apply(x, 2, function(arr) median(diff(sort(arr)))))
      mid.range <- mean(apply(x, 2, function(arr) diff(range(arr))))
      range.bd <- sort(c(5 * mid.diff, mid.range / 3))
      if(search == 'grid'){
        out <- expand.grid(
          poly.deg = seq(3),
          # kernel.bandwidth = 2 ^ (seq(len) - 3) * mid.diff * 10
          kernel.bandwidth = seq(range.bd[1], range.bd[2], length.out = len)
        )
      }else{
        out <- data.frame(
          poly.deg = sample.int(3, size = len, replace = TRUE),
          kernel.bandwidth = runif(len, range.bd[1], range.bd[2])

        )
      }

      return(out)
    }
    , fit = function(x, y, wts, param, lev, last, weights, classProbs){
      # a function that fits the model
      df <- as.data.frame(x)
      names(df) <- names(x)
      kernel.f <-
        getKernelFunc(kernel = 'gauss', bandwidth = param$kernel.bandwidth)
      loc.poly <- locPolyDefine(
        terms = names(x), degree = param$poly.deg, data = df,
        optns = list(kernel.f = kernel.f)
      )
      # sometimes there will be x & y in the x
      loc.poly$data <- list(x = df, y = factor(y))
      return(loc.poly)
    }
    , predict = function(modelFit, newdata, preProc = NULL, submodels = NULL){
      # the function that creates predictions
      w.df <- newdata[modelFit$terms]
      arr.pred.class <- rep(NA, nrow(w.df))
      # ls.glm <- list()
      for(i in seq(nrow(w.df))){
        tgt.pnt <- as.numeric(w.df[i, ])
        wt <- modelFit$weight.f(tgt.pnt)
        # only keep those with wt > .Machine$double.eps^(1/3)
        wt[wt < .Machine$double.eps^(1/3)] <- 0
        idx.keep <- wt != 0
        res.glm <- multiNomialGLM(
          y = modelFit$data$y[idx.keep],
          x = modelFit$mdl.mat.f(tgt.pnt, wt),
          wt = wt[idx.keep]
        )
        # fitted intercepts
        fitted.intercept <- res.glm$coefficients[1, ]
        lvl.nm <- names(fitted.intercept) # label for levels
        # recover class probability
        res.prob <- exp(fitted.intercept)
        res.prob <- res.prob / sum(res.prob)
        arr.pred.class[i] <- lvl.nm[which.max(res.prob)]
      }
      return(arr.pred.class)
    }
    , prob = function(modelFit, newdata, preProc = NULL, submodels = NULL){
      # browser()
      # a function that can be used to create class probabilities (if applicable)
      w.df <- newdata[modelFit$terms]
      n.lvl <- length(unique(modelFit$data$y))
      mat.prob <- matrix(NA, nrow = nrow(w.df), ncol = n.lvl)
      # ls.glm <- list()
      for(i in seq(nrow(w.df))){
        tgt.pnt <- as.numeric(w.df[i, ])
        wt <- modelFit$weight.f(tgt.pnt)
        # only keep those with wt > .Machine$double.eps^(1/3)
        wt[wt < .Machine$double.eps^(1/3)] <- 0
        idx.keep <- wt != 0
        res.glm <- multiNomialGLM(
          y = modelFit$data$y[idx.keep],
          x = modelFit$mdl.mat.f(tgt.pnt, wt),
          wt = wt[idx.keep]
        )
        # fitted intercepts
        fitted.intercept <- res.glm$coefficients[1, ]
        lvl.nm <- names(fitted.intercept) # label for levels
        # recover class probability
        res.prob <- exp(fitted.intercept)
        res.prob <- res.prob / sum(res.prob)
        mat.prob[i, ] <- res.prob
      }
      colnames(mat.prob) <- lvl.nm
      return(mat.prob)
    }
    , sort = function(x) x[order(x$poly.deg), ]
  )

  return(res)
}
