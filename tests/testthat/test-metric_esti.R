# Core functions for metric estimation

test_that('fitMetric and estiMetric', {

  d <- 3
  manifold <- spaceEuclidean(d)
  set.seed(1)
  obsv.coord <- manifold$genPnt(10^4)
  all.edge <- allEdge(obsv.coord, local.reach = 0.2)
  idx.edge <- all.edge[all.edge[, 1] != all.edge[, 2], , drop = F]
  idx.edge <- idx.edge[sample.int(nrow(idx.edge), 10^5), , drop = F]
  sqdist <- manifold$dist(
    obsv.coord[idx.edge[, 1], ], obsv.coord[idx.edge[, 2], ]
  ) ^ 2

  # input coord in data
  data.w.coord <- cbind(
    sqdist, obsv.coord[idx.edge[, 1], ], obsv.coord[idx.edge[, 2], ]
  )
  data.w.coord <- as.data.frame(data.w.coord)
  names(data.w.coord) <- c('y', sprintf('start%s', seq(d)), sprintf('end%s', seq(d)))
  formula.w.coord <- as.formula(sprintf(
    "y ~ (%s) : (%s)",
    paste(sprintf('start%s', seq(d)), collapse = ' + '),
    paste(sprintf('end%s', seq(d)), collapse = ' + ')
  ))

  # or use indices
  data.w.idx <- cbind(sqdist, idx.edge)
  data.w.idx <- as.data.frame(data.w.idx)
  names(data.w.idx) <- c('y', sprintf('p%s', seq(2)))
  formula.w.idx <- y ~ p1 : p2

  fit.w.coord <- fitMetric(formula.w.coord, data.w.coord)
  fit.w.idx <- fitMetric(formula.w.idx, data.w.idx, coord = obsv.coord)

  expect_equal(
    estiMetric(rep(0.5, d), fit.w.coord), estiMetric(rep(0.5, d), fit.w.idx)
  )

  set.seed(10)
  target <- matrix(runif(d * 10), ncol = d)
  true.metric <- apply(target, 1, manifold$metric, simplify = F)

  # interchangeable model
  expect_equal(
    estiMetric(target, fit.w.coord), estiMetric(target, fit.w.idx)
  )

  esti.metric <- estiMetric(target, fit.w.coord)
  names(esti.metric) <- NULL

  # check symmetric and positive definite
  expect_true(all(sapply(esti.metric, isSymmetric)))
  expect_true(all(sapply(esti.metric, function(metric){
    eig.val <- eigen(metric, symmetric = T, only.values = T)
    all(eig.val$values >= 0)
  })))
  # check estimation accuracy
  expect_equal(esti.metric, true.metric)

  ### binary censoring, not good for the moment, later TBDTBD
  dat.binary <- data.w.coord
  dat.binary$y <- sigmoid.f(20 * (dat.binary$y - 0.05))
  # hist(dat.binary$y)
  set.seed(1)
  dat.binary$y <- rbinom(nrow(dat.binary), 1, prob = dat.binary$y)
  # dat.binary$y %>% table()
  fit.binary <- fitMetric(formula.w.coord, dat.binary,optns = list(
    local.reach = 0.5, n.local = Inf, method.trim = 'proximity'
    , type = 'thre', intercept = TRUE
  ))
  set.seed(10)
  target <- manifold$genPnt(10)
  true.metric <- apply(target, 1, manifold$metric, simplify = F)
  true.metric <- lapply(true.metric, function(x) 20 * x)
  esti.metric <- estiMetric(target, fit.binary)
  names(esti.metric) <- NULL
  # check symmetric and positive definite
  expect_true(all(sapply(esti.metric, isSymmetric)))
  expect_true(all(sapply(esti.metric, function(metric){
    eig.val <- eigen(metric, symmetric = T, only.values = T)
    all(eig.val$values >= 0)
  })))
  # check estimation accuracy, not so well, honestly.
  expect_equal(esti.metric, true.metric, tolerance = 0.125)

  ### comparing edges
  set.seed(1)
  comp.edge <- genCompareGraph(c(500, 500), all.edge)
  sq.dist <- list(
    manifold$dist(
      obsv.coord[comp.edge[, 1], ], obsv.coord[comp.edge[, 2], ]
    ) ^ 2
    ,
    manifold$dist(
      obsv.coord[comp.edge[, 3], ], obsv.coord[comp.edge[, 4], ]
    ) ^ 2
  )
  arr.prob <- sigmoid.f(25 * do.call(`-`, sq.dist))
  # arr.prob %>% hist
  set.seed(1)
  dat.compare <- data.frame(
    y = rbinom(length(arr.prob), 1, prob = arr.prob),
    comp.edge
  )
  names(dat.compare) <- c('y', sprintf('p%s', seq(4)))
  # dat.compare %>% head
  fit.compare <- fitMetric(
    y ~ p1:p2:p3:p4, dat.compare, obsv.coord, optns = list(
      local.reach = 0.5
      , n.local = 10^4, method.trim = 'proximity'
      # , type = 'comp', intercept = FALSE
    )
  )

  # test guess model
  expect_identical(fit.compare$optns$type, 'compare')

  # expect_equal(
  #   estiMetric(rep(0.1, d), fit.compare)[[1]],
  #   diag(rep(25, d))
  #   , tolerance = 0.075
  # )

  set.seed(10)
  target <- manifold$genPnt(10)
  true.metric <- apply(target, 1, manifold$metric, simplify = F)
  true.metric <- lapply(true.metric, function(x) 25 * x)
  esti.metric <- estiMetric(target, fit.compare)
  names(esti.metric) <- NULL
  # check symmetric and positive definite
  expect_true(all(sapply(esti.metric, isSymmetric)))
  expect_true(all(sapply(esti.metric, function(metric){
    eig.val <- eigen(metric, symmetric = T, only.values = T)
    all(eig.val$values >= 0)
  })))
  # check estimation accuracy, not so well, honestly.
  expect_equal(esti.metric, true.metric, tolerance = 0.05)

})

test_that("locMetric", {

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
    locMetric(arr.dist ^ 2, coord.diff)$mat.g
  )
  # normal noise (after squared)
  expect_equal(
    diag(rep(1, d)),
    locMetric(
      arr.dist ^ 2 + rnorm(length(arr.dist), sd = sd(arr.dist) / 50),
      coord.diff
    )$mat.g
    , tolerance = 5e-3
  )
  # binary thresholding
  expect_equal(
    diag(rep(1, d)),
    locMetric(
      arr.thre, coord.diff
      , optns = list(type = 't', intercept = FALSE)
    )$mat.g
    , tolerance = 0.065 #0.275
  )
  # binary thresholding w/ intercept
  expect_equal(
    diag(rep(1, d)),
    locMetric(
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
    locMetric(
      arr.comp, ls.coord.diff
      , optns = list(type = 'comp', intercept = FALSE)
    )$mat.g
    , tolerance = 0.035 #0.175
  )
  # comparing distance w/ intercept
  res.itcpt <- locMetric(
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
