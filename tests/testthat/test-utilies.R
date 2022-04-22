# Utility functions

test_that("QuadForm", {
  arr.x <- seq(3)
  expect_identical(
    c(1, 4, 9, 2, 3, 6),
    as.numeric(vecQuad(arr.x))
  )

  mat.x <- matrix(runif(3 * 4), ncol = 3)
  idx.col <- vecQuadIdx(3)
  expect_identical(
    matrix(c(
      rep(seq(3L), each = 2),
      c(1L, 2L, 1L, 3L, 2L, 3L)
    ), nrow = 2)
    , idx.col
  )

  man.quad.vec <-
    mat.x[, idx.col[1, ]] * mat.x[, idx.col[2, ]]
  fun.quad.vec <- vecQuad(mat.x)

  expect_identical(
    man.quad.vec,
    matrix(fun.quad.vec, nrow = nrow(mat.x))
  )

  d <- 4L
  mat.x <- matrix(runif(d * 10), ncol = d)
  df.x <- as.data.frame(mat.x)
  names(df.x) <- sprintf('x%s', seq(d))
  mat.a <- matrix(runif(d*d), ncol = d)
  mat.a <- (mat.a + t(mat.a)) / 2
  vec.a <-
    c(diag(mat.a), 2 * mat.a[lower.tri(mat.a)])
  # dimension check
  expect_identical(
    d * (d + 1) / 2,
    as.numeric(ncol(vecQuad(mat.x)))
  )
  expect_identical(
    rep(nrow(mat.x), 2),
    dim(quadForm(t(mat.x), mat.a))
  )
  # quadratic form check
  expect_equal(
    as.numeric(vecQuad(mat.x) %*% vec.a),
    diag(quadForm(t(mat.x), mat.a))
  )

})


test_that("Sigmoid and logit function", {
  arr.x <- seq(-15, 15, 0.01)
  arr.p <- seq(1e-5, 1-(1e-5), 0.01)
  rbuildin <- stats::binomial()
  expect_equal(
    arr.x,
    sigmoid.f.inv(sigmoid.f(arr.x))
  )
  expect_equal(
    arr.p,
    sigmoid.f(sigmoid.f.inv(arr.p))
  )
  expect_equal(
    rbuildin$linkinv(arr.x),
    sigmoid.f(arr.x)
  )
  expect_equal(
    rbuildin$linkfun(arr.p),
    sigmoid.f.inv(arr.p)
  )
})

test_that("Square neighborhood identification", {
  set.seed(1)
  d <- 3
  coord <- matrix(runif(d * 2000), ncol = d)
  edge.res <- allEdge(coord, 0.1)
  locW.res <- locWindow(coord, coord, 0.1, 'matrix')
  expect_identical(
    locW.res[locW.res[, 1] <= locW.res[, 2], ],
    edge.res
  )

  # microbenchmark::microbenchmark(list = alist(
  #   locW = {
  #     locW.res <- locWindow(coord, coord, 0.1, 'matrix')
  #     locW.res[locW.res[, 1] <= locW.res[, 2], ]
  #   },
  #   edge = allEdge(coord, 0.1)
  # ))

})


test_that("edge data format conversion", {
  set.seed(1)
  d <- 3
  coord <- matrix(runif(d * 10), ncol = d)
  edge.mat <- locWindow(rep(0, d), coord, 0.5, 'matrix')
  edge.list <- locWindow(rep(0, d), coord, 0.5, 'list')
  expect_identical(edgeMat2List(edge.mat), edge.list)
  expect_identical(edgeList2Mat(edge.list), edge.mat)
  expect_identical(edgeMat2List(edgeList2Mat(edge.list)), edge.list)
  expect_identical(edgeList2Mat(edgeMat2List(edge.mat)), edge.mat)
  expect_identical(formatGraph(edge.mat), edge.list)
  expect_identical(formatGraph(edge.list), edge.mat)
  expect_identical(formatGraph(edge.mat, 'matrix'), edge.mat)
  expect_identical(formatGraph(edge.list, 'list'), edge.list)
})
