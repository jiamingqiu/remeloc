test_that("Tensor vectorization with symmetries", {

  d <- 3
  t.metric <- termMetric(d = d)
  t.metric$df.idx
  tst <- matrix(rnorm(d^2), d, d)
  tst <- tst + t(tst)
  expect_identical(
    with(t.metric, tst %>% tsr2vec() %>% vec2tsr),
    tst
  )

  t.pchrist <- termPreChrist(d = d)
  t.pchrist$df.idx
  tst <- array(rnorm(d^3), rep(d, 3))
  for(k in seq(d)){ tst[,, k] <- tst[,, k] + t(tst[,, k]) }
  expect_identical(
    with(t.pchrist, tst %>% tsr2vec() %>% vec2tsr),
    tst
  )
  # tst %>% t.pchrist$tsr2vec() %>% length
  # tst %>% length

})

test_that("Tensor contraction", {

  set.seed(1)
  # tensor action
  d <- c(2, 3, 4, 5)
  tst.04 <- array(runif(prod(d)), dim = d)
  ls.vec <- list()
  for(i in seq_along(d)){
    ls.vec[[i]] <- runif(d[i])
  }
  act.pck <- do.call(getTensorAction, c(list(tst.04), ls.vec))

  act.for <- tst.04
  for(i in seq(d[1])){ for(j in seq(d[2])){
    for(k in seq(d[3])){ for(l in seq(d[4])){
      act.for[i, j, k, l] <- act.for[i, j, k, l] *
        ls.vec[[1]][i] * ls.vec[[2]][j] * ls.vec[[3]][k] * ls.vec[[4]][l]
    }}
  }}
  act.for <- sum(act.for)
  expect_identical(act.for, act.pck)

  # sectional curvature
  d <- 2
  manifold <- spaceHyperbolic(d = d, r = 5, model = 'ball')
  cvt04.f <- get04Curvature(manifold$metric, d = d)
  cvt.act <- getTensorAction(cvt04.f)
  metric.act <- getTensorAction(manifold$metric)
  # target <- rep(1, d)
  sec.cvt <- function(target) {
    cvt.act(target, c(1, 0), c(0, 1), c(0, 1), c(1, 0)) / (
      metric.act(target, c(1, 0), c(1, 0)) *
        metric.act(target, c(0, 1), c(0, 1)) -
        metric.act(target, c(1, 0), c(0, 1)) ^ 2
    )
  }
  expect_equal(
    apply(manifold$genPnt(10), 1, sec.cvt),
    rep(-1 / 5^2, 10), tolerance = 5e-6
  )

  # sectional curvature
  sec.cvt.f <- getSecCurvature(manifold$metric, d = d)
  expect_equal(
    apply(manifold$genPnt(10), 1, sec.cvt.f, v = c(1, 0), w = c(0, 1)),
    rep(-1 / 5^2, 10), tolerance = 5e-6
  )
  set.seed(1)
  expect_equal(
    apply(manifold$genPnt(10), 1, sec.cvt.f, v = c(1, 0), w = c(0, 1)),
    apply(manifold$genPnt(10), 1, sec.cvt)
    , tolerance = 5e-6
  )

  # higher dimension sphere
  d <- 5
  manifold <- spaceSphere(d = d, r = 5)
  sec.cvt.f <- getSecCurvature(manifold$metric, d = d)
  expect_equal(
    apply(
      manifold$genPnt(10), 1, sec.cvt.f,
      v = c(1, rep(0, d - 1)), w = c(rep(0, d - 1), 1)
    ),
    rep(1 / 5^2, 10), tolerance = 3e-3
  )


})

test_that("vectorizing symmetric matrix", {

  mat <- matrix(runif(3 * 3), 3, 3)
  mat <- (mat + t(mat)) / 2
  vec <- c(diag(mat), mat[lower.tri(mat)])

  expect_identical(vec2SymMat(vec), mat)
  expect_identical(symMat2Vec(mat), vec)
  expect_identical(mat[vecQuadIdx(3, T)], vec)

})



