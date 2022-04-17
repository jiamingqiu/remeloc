# Core functions for metric estimation

test_that("locMetric.g", {

  d <- 3

  ## Euclidean space
  set.seed(1)
  # coord of points
  coord.pnts <- matrix(runif(1000 * d, 0, 2), ncol = d)
  # pairing points
  idx.dist <- t(utils::combn(nrow(coord.pnts), 2))
  idx.dist <- idx.dist[
    sample(nrow(idx.dist), size = 25000),
  ]
  coord.diff <-
    coord.pnts[idx.dist[, 1], ] - coord.pnts[idx.dist[, 2], ]
  arr.dist <- sqrt(rowSums(coord.diff ^ 2))
  prob.thre <- sigmoid.f(arr.dist ^ 2)
  prob.thre.itcpt <- sigmoid.f(arr.dist ^ 2 - mean(arr.dist ^ 2))
  arr.thre <- rbinom(length(arr.dist), 1, prob = prob.thre)
  arr.thre.itcpt <- rbinom(length(arr.dist), 1, prob = prob.thre.itcpt)

  # # rmk: for thresholding and compare type, check the prob so cover most of 0-1
  # # otherwise esti. will be very poor despite num_obsv.
  # hist(arr.dist)
  # hist(prob.thre)
  # hist(prob.thre.itcpt)

  # noiseless observation
  expect_equal(
    diag(rep(1, d)),
    locMetric.g(arr.dist ^ 2, coord.diff)$mat.g
  )
  # normal noise (after squared)
  expect_equal(
    diag(rep(1, d)),
    locMetric.g(
      arr.dist ^ 2 + rnorm(length(arr.dist), sd = sd(arr.dist) / 50),
      coord.diff
    )$mat.g
    , tolerance = 5e-3
  )
  # binary thresholding
  expect_equal(
    diag(rep(1, d)),
    locMetric.g(
      arr.thre, coord.diff
      , optns = list(type = 't', intercept = FALSE)
    )$mat.g
    , tolerance = 0.065 #0.275
  )
  # binary thresholding w/ intercept
  expect_equal(
    diag(rep(1, d)),
    locMetric.g(
      arr.thre.itcpt, coord.diff
      , optns = list(type = 't', intercept = TRUE)
    )$mat.g
    , tolerance = 0.045
  )

  # comparing distance
  set.seed(1)
  idx.comp <- t(utils::combn(
    sample(length(arr.dist), 1000), 2)
  )
  idx.comp <- idx.comp[
    sample(nrow(idx.comp), size = 10000),
  ]
  ls.coord.diff <- list(
    coord.diff[idx.comp[, 1], ],
    coord.diff[idx.comp[, 2], ]
  )
  prob.comp <- sigmoid.f(
    arr.dist[idx.comp[, 1]]^2 - arr.dist[idx.comp[, 2]]^2
  )
  arr.comp <- rbinom(length(prob.comp), 1, prob = prob.comp)

  # hist(prob.comp)
  # table(arr.comp)

  # comparing distance
  expect_equal(
    diag(rep(1, d)),
    locMetric.g(
      arr.comp, ls.coord.diff
      , optns = list(type = 'comp', intercept = FALSE)
    )$mat.g
    , tolerance = 0.035 #0.175
  )
  # comparing distance w/ intercept
  res.itcpt <- locMetric.g(
    arr.comp, ls.coord.diff
    , optns = list(type = 'comp', intercept = TRUE)
  )
  expect_equal(
    diag(rep(1, d)),
    res.itcpt$mat.g
    , tolerance = 0.035 #0.175
  )
  # the fitted intercept should be close to 0
  expect_equal(
    0, as.numeric(res.itcpt$loc.fit$coefficients[1])
    , tolerance = 0.025 #0.007
  )

  # expect_true(FALSE)
})
