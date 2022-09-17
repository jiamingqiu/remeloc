
mat2PD <- function(mat){
  # truncate a symmetric matrix to be non-negative definite
  stopifnot(isSymmetric(mat))
  eig <- eigen(mat, symmetric = T)
  if(all(eig$values >= .Machine$double.eps^(1/3)))
    return(mat)

  eig$values[eig$values < .Machine$double.eps^(1/3)] <-
    .Machine$double.eps^(1/3)

  return(
    tcrossprod(
      tcrossprod(eig$vectors, diag(eig$values)),
      eig$vectors
    )
  )
}



##### Model terms related ######################################################

termMetric <- function(d, method = 'adhoc'){
  # return tools(functions) used to construct metric term
  # sum_{i, j} (cd0[i] - cd1[i]) * (cd0[j] - cd1[j]) * metric[i, j]

  method <- match.arg(method, c('adhoc', 'basis'))

  if(method == 'adhoc'){

    # index for vectorizing curvature tensor action
    df.idx <- data.frame(idx = vecQuadIdx(d, idx.vec = T))
    df.idx <- cbind(df.idx, as.data.frame(t(vecQuadIdx(d))))
    names(df.idx) <- c('idx', sprintf('i%s', seq(2)))
    df.idx$off.diag <- T
    df.idx$off.diag[seq(d)] <- F
    df.idx$subscript <- base::Reduce(
      function(x, y) sprintf('%s,%s', x, y),
      x = df.idx[, sprintf('i%s', seq(2)), drop = F]
    )

    # # basis for the vectorization
    # vec.tsr.basis <- matrix(0, nrow = nrow(df.idx), ncol = d * (d + 1) / 2)
    # vec.tsr.basis[df.idx$idx, df.idx$idx] <- 1
    vec.tsr.basis <- NULL

    # model matrix constructor
    modelMat <- function(cd.cd, cd0, cd1) {
      if(missing(cd.cd)) cd.cd <- cd0 - cd1
      mat <- cd.cd[, df.idx$i1, drop = F] * cd.cd[, df.idx$i2, drop = F]
      mat[, -seq(d)] <- 2 * mat[, -seq(d), drop = F]
      colnames(mat) <- sprintf('metric.%s', seq(ncol(mat)))
      return(mat)
    }
    # tensor <--> vector translators
    tsr2vec <- function(tensor) {
      # compute coefficient from input array
      # symMat2Vec(tensor)
      return(tensor[df.idx$idx])
    }
    vec2tsr <- function(vec){
      return(vec2SymMat(vec))
    }

  }else{

    # index for vectorizing curvature tensor action
    df.idx <- vecTensorIdx(d = d, k = 2, sym.eq = 'ij=ji')
    # basis for the vectorization
    vec.tsr.basis <- as.matrix(df.idx[
      , stringr::str_detect(names(df.idx), 'tsr\\.e\\d+'), drop = F
    ])
    # model matrix constructor
    modelMat <- function(cd.cd, cd0, cd1) {
      mat <- modelMatMetric(df.idx, cd.cd, cd0, cd1)
      colnames(mat) <- sprintf('metric.%s', seq(ncol(mat)))
      return(mat)
    }

    # tensor <--> vector translators
    tsr2vec <- function(tensor) {
      # compute projection score from input array, not vectorized!
      colSums(tensor[df.idx$idx][row(vec.tsr.basis)] * vec.tsr.basis)
    }
    vec2tsr <- function(vec){
      tsr <- array(0, dim = rep(d, 2))
      tsr[df.idx$idx] <- rowSums(vec.tsr.basis * vec[col(vec.tsr.basis)])
      return(tsr)
    }

  }
  return(list(
    df.idx = df.idx, vec.tsr.basis = vec.tsr.basis,
    modelMat = modelMat, tsr2vec = tsr2vec, vec2tsr = vec2tsr
  ))

}

termPreChrist <- function(d, method = 'adhoc'){
  # return tools(functions) used to construct pre-Christ term
  # sum_{i, m, n} ( cd0[i] - cd1[i] ) * (cd0[m] * cd0[n] - cd1[m] * cd1[n]) * (
  #   Christ_{m, n}^j * metric[i, j]
  # )

  method <- match.arg(method, c('adhoc', 'basis'))

  if(method == 'adhoc'){

    len.2.tsr <- d ^ 2
    # len.2.tsr.vec <- d * (d + 1) / 2

    # index for vectorizing curvature tensor action
    df.idx <- data.frame(idx = vecQuadIdx(d, idx.vec = T))
    df.idx <- cbind(df.idx, as.data.frame(t(vecQuadIdx(d))))
    names(df.idx) <- c('idx', sprintf('i%s', seq(2)))
    df.idx <- lapply(seq(d), function(idx.3rd){
      tm <- df.idx
      tm$idx <- tm$idx + (idx.3rd - 1) * len.2.tsr
      tm$i3 <- idx.3rd
      tm$off.diag <- T
      tm$off.diag[seq(d)] <- F
      tm
    })
    df.idx <- do.call(rbind, df.idx)
    df.idx$subscript <- base::Reduce(
      function(x, y) sprintf('%s,%s', x, y),
      x = df.idx[, sprintf('i%s', seq(3))]
    )

    # browser();QWER

    # # basis for the vectorization
    # vec.tsr.basis <- matrix(0, nrow = nrow(df.idx), ncol = d * (d + 1) / 2)
    # vec.tsr.basis[df.idx$idx, df.idx$idx] <- 1
    vec.tsr.basis <- NULL

    # model matrix constructor
    modelMat <- function(cd.cd, cd0, cd1) {
      if(missing(cd.cd)) cd.cd <- cd0 - cd1
      mat <- with(df.idx, {
        cd.cd[, i3, drop = F] * (
          cd0[, i1, drop = F] * cd0[, i2, drop = F] -
            cd1[, i1, drop = F] * cd1[, i2, drop = F]
        )
      })
      mat[, df.idx$off.diag] <- 2 * mat[, df.idx$off.diag, drop = F]
      colnames(mat) <- sprintf('pre.christ.%s', seq(ncol(mat)))
      return(mat)
    }
    # tensor <--> vector translators
    tsr2vec <- function(tensor) {
      # compute coefficient from input array
      return(tensor[df.idx$idx])
    }
    vec2tsr <- function(vec){
      # browser();QWER
      # split into d 2-tensors
      ls.2.tsr <- split(vec, df.idx$i3)
      # each into symmetric matrix
      res <- lapply(ls.2.tsr, vec2SymMat)
      return(array(do.call(c, res), dim = rep(d, 3)))
    }

  }else{

    # index for vectorizing curvature tensor action
    df.idx <- vecTensorIdx(d = d, k = 3, sym.eq = 'ijk=jik')
    # basis for the vectorization
    vec.tsr.basis <- as.matrix(df.idx[
      , stringr::str_detect(names(df.idx), 'tsr\\.e\\d+'), drop = F
    ])
    # model matrix constructor
    modelMat <- function(cd.cd, cd0, cd1) {
      mat <- modelMatPreChrist(df.idx, cd.cd, cd0, cd1)
      colnames(mat) <- sprintf('pre.christ.%s', seq(ncol(mat)))
      return(mat)
    }

    # tensor <--> vector translators
    tsr2vec <- function(tensor) {
      # compute projection score from input array, not vectorized!
      colSums(tensor[df.idx$idx][row(vec.tsr.basis)] * vec.tsr.basis)
    }
    vec2tsr <- function(vec){
      tsr <- array(0, dim = rep(d, 3))
      tsr[df.idx$idx] <- rowSums(vec.tsr.basis * vec[col(vec.tsr.basis)])
      return(tsr)
    }

  }

  return(list(
    df.idx = df.idx, vec.tsr.basis = vec.tsr.basis,
    modelMat = modelMat, tsr2vec = tsr2vec, vec2tsr = vec2tsr
  ))

}

#' Wrapper for pre-Christ to Christoffel symbol
#'
#' @param pre.christ a (0, 3)-tensor
#' @param metric the metric tensor
#'
#' @return either a (1, 2)-tensor or a function (if inputs are function)
#' computing it.
#' @export
#'
#' @details The pre-Christ is defined as
#' \deqn{F_{mni} = \Gamma_{mn}^j g_{ij},}
#' where \eqn{\Gamma_{mn}^j} is the Christoffel symbol, and \eqn{g_{ij}} is the
#' metric tensor. The return follows \code{F[m, n, i] = }\eqn{F_{mni}}.
preChrist2Christ <- function(pre.christ, metric){

  # computes Christoffel symbol from pre.christ

  if(inherits(pre.christ, 'function')){
    stopifnot(inherits(metric, 'function'))
    res.f <- function(x) preChrist2Christ(pre.christ(x), metric(x))
    return(res.f)
  }

  d <- nrow(metric)

  # next, arr.christ[m,n,i,j] = preChrist[m, n, i] = Christ[m,n,k] * metric[i,k]
  arr.christ <- array(pre.christ, dim = rep(d, 4))
  # next, arr.christ[m,n,i,j] = preChrist[m, n, i] * metric_inv[i, j]
  arr.christ <- arr.christ * remeloc::reindex(
    array(MASS::ginv(metric), dim = rep(d, 4)),
    from = 'ijmn', to = 'mnij' # reidx[m,n,i,j] = old[i,j,m,n] = metric_inv[i,j]
  )
  # next contraction for sum_i(arr.christ[m,n,i,j])
  arr.christ <- rowSums(
    remeloc::reindex(arr.christ, from = 'mnij', to = 'mnji'),
    dims = 3
  )

  return(arr.christ)

}

termRiemCvt <- function(d){
  # return tools(functions) used to construct pre-Christ term
  # sum_{i, m, n} ( cd0[i] - cd1[i] ) * (cd0[m] * cd0[n] - cd1[m] * cd1[n]) * (
  #   Christ_{m, n}^j * metric[i, j]
  # )

  # index for vectorizing curvature tensor action
  df.idx <- vecRiemCvtIdx(d = d)
  # basis for the vectorization
  vec.tsr.basis <- as.matrix(df.idx[
    , stringr::str_detect(names(df.idx), 'tsr\\.e\\d+'), drop = F
  ])
  # model matrix constructor
  modelMat <- function(cd.cd, cd0, cd1) {
    mat <- -1/3 * modelMatRiemCvt(df.idx, cd0, cd1, cd1, cd0)
    colnames(mat) <- sprintf('cvt.%s', seq(ncol(mat)))
    return(mat)
  }

  # tensor <--> vector translators
  tsr2vec <- function(tensor) {
    # compute projection score from input array
    colSums(tensor[df.idx$idx][row(vec.tsr.basis)] * vec.tsr.basis)
  }
  vec2tsr <- function(vec){
    tsr <- array(0, dim = rep(d, 4))
    tsr[df.idx$idx] <- rowSums(vec.tsr.basis * vec[col(vec.tsr.basis)])
    return(tsr)
  }

  return(list(
    df.idx = df.idx, vec.tsr.basis = vec.tsr.basis,
    modelMat = modelMat, tsr2vec = tsr2vec, vec2tsr = vec2tsr
  ))

}

modelMatMetric <- function(df.idx, cd.cd, cd0, cd1){

  # computes model matrix of the metric term
  # cd0 and cd1: coord diff of either endpoints to target point
  # cd.cd: cd0 - cd1

  # ls.cd <- list(...)
  # stopifnot(length(ls.cd) <= 2)
  # if(length(ls.cd) == 2)
  #   cd.cd <- ls.cd[[1]] - ls.cd[[2]]
  # else
  #   cd.cd <- ls.cd[[1]]

  if(missing(cd.cd))
    cd.cd <- cd0 - cd1

  arr.r <- with(df.idx, {
    cd.cd[, i1] * cd.cd[, i2]
  }) # arr.r[, j] = term in df.idx[j, ]

  vec.tsr.basis <- as.matrix(
    df.idx[, stringr::str_detect(names(df.idx), 'tsr\\.e\\d+'), drop = F]
  )
  # computes arr.r %*% vec.tsr.basis, so the model
  # arr.r %*% vec.tsr = arr.r %*% vec.tsr %*% proj.score
  # mat.r <- apply(vec.tsr.basis, 2, function(x) {
  #   rowSums(arr.r * x[col(arr.r)])
  # })
  mat.r <- tcrossprod(arr.r, t(vec.tsr.basis))

  return(mat.r)
}

modelMatPreChrist <- function(df.idx, cd.cd, cd0, cd1){

  # computes model matrix of the pre-Christ term, i.e., the
  # sum_{i, m, n} ( cd0[i] - cd1[i] ) * (cd0[m] * cd0[n] - cd1[m] * cd1[n]) * (
  #   Christ_{m, n}^j * metric[i, j]
  # )
  # df.idx: index data frame for symmetries
  # ...: 2 arguements, cd0 and cd1, the difference in coordinates from either
  # endpoints to the target.

  # ls.cd <- list(...)
  # stopifnot(length(ls.cd) <= 3)
  # cd.cd <- ls.cd[[1]] - ls.cd[[2]]

  if(missing(cd.cd))
    cd.cd <- cd0 - cd1

  arr.r <- with(df.idx, {
    (
      cd0[, i1] * cd0[, i2] - cd1[, i1] * cd1[, i2]
    ) * cd.cd[, i3]
  }) # arr.r[, j] = term in df.idx[j, ]

  vec.tsr.basis <- as.matrix(
    df.idx[, stringr::str_detect(names(df.idx), 'tsr\\.e\\d+'), drop = F]
  )
  # browser();QWER
  # all.equal(
  #   tcrossprod(arr.r, t(vec.tsr.basis)),
  #   apply(vec.tsr.basis, 2, function(x) {
  #     rowSums(arr.r * x[col(arr.r)])
  #   })
  # )
  # microbenchmark::microbenchmark(list = alist(
  #   tcrossprod = tcrossprod(arr.r, t(vec.tsr.basis)), # mutch faster
  #   apply = apply(vec.tsr.basis, 2, function(x) {
  #     rowSums(arr.r * x[col(arr.r)])
  #   })
  # ))
  # computes arr.r %*% vec.tsr.basis, so the model
  # arr.r %*% vec.tsr = arr.r %*% vec.tsr %*% proj.score
  # mat.r <- apply(vec.tsr.basis, 2, function(x) {
  #   rowSums(arr.r * x[col(arr.r)])
  # })
  mat.r <- tcrossprod(arr.r, t(vec.tsr.basis))

  return(mat.r)
}

modelMatRiemCvt <- function(df.idx, x, y, z, w){

  # construct design matrix for local regression on Riemann (0,4)-curvature
  # df.idx: result from vecRiemCvtIdx
  # ...: 4 tensor arrays, no sanity check, make sure dimension match!
  # also, no adjustment for potential bias due to using coordinate

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
    df.idx[, stringr::str_detect(names(df.idx), 'tsr\\.e\\d+'), drop = F]
  )
  # put the coef for basis into columns
  # mat.r <- -1 / 3 * apply(vec.cvt.basis, 2, function(x) {
  #   rowSums(arr.r * x[col(arr.r)])
  # })
  mat.r <- - 1 / 3 * apply(vec.cvt.basis, 2, function(x) {
    rowSums(arr.r * x[col(arr.r)])
  })

  return(mat.r)

}

modelMatTensor <- function(df.idx, ...){

  # computes model matrix of a class of tensor (from idx data.frame)
  # acting on ...
  # for example, compute model matrix of
  # sum_{i,j,k,l} R[i,j,k,l] * w[i] * x[j] * y[k] * z[l]
  # according to df.idx

  ls.arr <- list(...)
  stopifnot(
    length(ls.arr) == length(stringr::str_detect(names(df.idx), 'i\\d+'))
  )
  arr.r <- with(df.idx, eval(parse(text = paste(
    sprintf('ls.arr[[%s]][, i%s]', seq_along(ls.arr)), collapse = ' * '
  )))) # arr.r[, j] = term in df.idx[j, ]

  vec.tsr.basis <- as.matrix(
    df.idx[, stringr::str_detect(names(df.idx), 'tsr\\.e\\d+'), drop = F]
  ) # 1col = 1 basis

  # computes arr.r %*% vec.tsr.basis, so the model
  # arr.r %*% vec.tsr = arr.r %*% vec.tsr %*% proj.score
  # mat.r <- apply(vec.tsr.basis, 2, function(x) {
  #   rowSums(arr.r * x[col(arr.r)])
  # })
  mat.r <- tcrossprod(arr.r, t(vec.tsr.basis))

  return(mat.r)
}

##### vectorization of Riemann curvature tensor ################################

getTensorAction <- function(tensor, ...){
  # tensor action on vectors

  if(inherits(tensor, 'function')){
    tensor <- tensor
    res.f <- function(target, ...) getTensorAction(tensor(target), ...)
    return(res.f)
  }

  stopifnot(inherits(tensor, 'array'))
  d.tensor <- dim(tensor)
  # stopifnot(length(d.tensor) >= 2)
  ls.vec <- list(...)
  stopifnot(length(ls.vec) == length(d.tensor))
  stopifnot(all(sapply(ls.vec, length) == d.tensor))

  # compute sum_ijkl tensor[i, j, k, l, ...] * w[i] * x[j] * y[k] * z[l] * ...
  tsr <- Reduce(
    function(tsr, v)  colSums(tsr * v, dims = 1),
    ls.vec[seq(length(ls.vec) - 1)], init = tensor
  )
  # the last is taken out to avoid colSums error
  return(sum(tsr * ls.vec[[length(ls.vec)]]))

}

#' Index table for vectorizing algebraic curvature tensor
#'
#' @param d dimension of vector space
#'
#' @return a data frame
#' @export
#'
#' @details See \code{\link{vecTensorIdx}} for details.
#' See p212, Lee (2018) for definition.
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
#' @param drop.zero whether to drop subscripts of always zero elements
#'
#' @return a data frame
#' @export
#'
#' @details The each row in the returned data frame represents one element in
#' the desired tensor, and with the following columns:
#' \describe{
#'   \item{\code{idx}}{The index of the element if flatten.}
#'   \item{\code{i1, ..., ik}, and \code{subscript}}{Subscripts.}
#'   \item{\code{tsr.e1, ...}}{Basis used for representing.}
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
  }else if(is.null(sym.eq)) {
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
  names(df.basis) <- sprintf('tsr.e%s', seq(ncol(df.basis)))
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
  idx <- vecQuadIdx(d, idx.vec = T)
  # idx <- idx[2, ] + d * (idx[1, ] - 1)
  return(as.numeric(mat[idx]))
}
# symMat2Vec(diag(3))
# symMat2Vec(matrix(seq(9), 3) + t(matrix(seq(9), 3)))

vec2SymMat <- function(vec){
  # construct symmetric matrix following vecQuadIdx
  d <- (sqrt(8 * length(vec) + 1) - 1) / 2
  mat <- matrix(NA, d, d)
  mat[vecQuadIdx(d, idx.vec = T)] <- vec
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
