# Core functions for metric estimation

test_that('options and model guessing', {

  expect_setequal(
    names(guess_model()), c(
      'method.trim', 'n.local', 'tensor.only', 'type', 'intercept',
      'deriv'
    )
  )
  expect_identical(set_optns(method.trim = 'random')$method.trim, 'random')
  expect_identical(
    guess_model(formula = y ~ p1:p2)[c('type', 'intercept')],
    list(type = 'distance', intercept = FALSE)
  )
  expect_identical(
    guess_model(formula = y ~ (x1 + x2):(z1 + z2))[c('type', 'intercept')],
    list(type = 'distance', intercept = FALSE)
  )
  expect_identical(
    guess_model(formula = y ~ p1:p2:p3:p4)[c('type', 'intercept')],
    list(type = 'compare', intercept = FALSE)
  )
  expect_identical(
    guess_model(
      formula = y ~ (e1x1 + e1x2):(e2x1 + e2x2):(e3x1 + e3x2):(e4x1 + e4x2)
    )[c('type', 'intercept')],
    list(type = 'compare', intercept = FALSE)
  )
  expect_error(
    guess_model(
      formula = y ~ (e1x1 + e1x2):(e2x1 + e2x2):(e3x1 + e3x2)
    )
  )

})

test_that('Mathematical accuracy of estimated tensors', {

  # skip_on_cran()
  # rather slow.

  d <- 3
  manifold <- spaceHyperbolic(d = d, r = 1, model = 'ball')

  christ.f <- getChristoffel(manifold$metric, d = d)
  cvt04.f <- get04Curvature(manifold$metric, d = d)

  set.seed(1)
  obsv.coord <- manifold$genPnt(10^6)
  system.time({
    # polar/stereo/poincare ball
    graph <- genGraph(5 * 10^5, obsv.coord, local.reach = 0.05)
    # graph <- genGraph(5 * 10^5, obsv.coord, local.reach = 0.1) # stereographic
    # graph <- genGraph(5 * 10^5, obsv.coord, local.reach = 0.01) # halfplane
  }) # ~ 60s for 5 * 10^5
  # target <- manifold$genGrid(32) * 0.75 + 0.25 / 2 # halfplane
  # target <- manifold$genGrid(32) # stereographic
  target <- manifold$genGrid(16) * 0.75 # poincare ball
  # target <- manifold$genGrid(32) * 0.75 + (pi * (1 - 0.75)) / 2 # polar
  # plt.dat <- plotGraph(obsv.coord, graph, max.n.edge = 500) +
  #   geom_point(
  #     data = as.data.frame(target) %>% setNames(c('x', 'y')),
  #     aes(x, y), colour = 'blue', size = 2
  #   )
  # plt.dat + coord_fixed()
  dat <- data.frame(
    sqd = manifold$dist(obsv.coord[graph[, 1], ], obsv.coord[graph[, 2], ]) ^ 2
  )
  dat <- cbind(dat, as.data.frame(graph))
  names(dat) <- c('sqd', 'p1', 'p2')
  # set.seed(1)
  # dat$y <- pmax(0, dat$sqd + rnorm(nrow(dat)) * 0.1 * dat$sqd)
  dat$y <- dat$sqd

  fit <- fitMetric(y ~ p1 : p2, data = dat, coord = obsv.coord, optns = list(
    local.reach = 0.2, deriv = 2, tensor.only = T
  ))
  # system.time({
    ls.esti <- estiMetric(target, fit)
  # }) # ~ 55s

  # compute relative Frobenius error
  relFrobErr <- function(true, esti){
    sum( (true - esti) ^ 2 ) / sum( true ^ 2)
  }
  expect_true(all(mapply(
    relFrobErr,
    true = apply(target, 1, manifold$metric, simplify = F),
    esti = ls.esti$metric
  ) <= 0.004))
  expect_true(all(mapply(
    relFrobErr,
    true = apply(target, 1, getChristoffel(manifold$metric, d), simplify = F),
    esti = ls.esti$christ
  ) <= 0.05))

  # # well, you can see there must be additional bias in the curvature term. TBD!
  # mapply(
  #   relFrobErr,
  #   true = apply(target, 1, cvt04.f, simplify = F),
  #   esti = ls.esti$curvature
  # ) %>% hist


})

test_that('fitMetric and estiMetric', {

  d <- 3
  manifold <- spaceEuclidean(d)
  # manifold <- spaceHyperbolic(d, r = 5)
  set.seed(1)
  obsv.coord <- manifold$genPnt(10^4)
  all.edge <- allEdge(obsv.coord, local.reach = 0.2)
  idx.edge <- all.edge[all.edge[, 1] != all.edge[, 2], , drop = F]
  idx.edge <- idx.edge[sample.int(nrow(idx.edge), 5 * 10^5), , drop = F]
  sqdist <- manifold$dist(
    obsv.coord[idx.edge[, 1], ], obsv.coord[idx.edge[, 2], ]
  ) ^ 2
  # sqdist %>% hist
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
  # target <- matrix(runif(d * 10), ncol = d)
  target <- manifold$genPnt(10) * 0.8
  true.metric <- apply(target, 1, manifold$metric, simplify = F)

  # interchangeable model
  expect_equal(
    estiMetric(target, fit.w.coord), estiMetric(target, fit.w.idx)
  )

  fit <- fitMetric(formula.w.idx, data.w.idx, coord = obsv.coord, optns = list(
    local.reach = 0.2, n.local = c(10, 10^5)
  ))
  esti.metric <- estiMetric(target, fit)
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
    local.reach = 0.2, n.local = Inf, method.trim = 'proximity'
    , type = 'thre', intercept = TRUE
  ))
  true.metric.binary <- lapply(true.metric, function(x) 20 * x)
  esti.metric <- estiMetric(target, fit.binary)$metric
  names(esti.metric) <- NULL
  # check symmetric and positive definite
  expect_true(all(sapply(esti.metric, isSymmetric)))
  expect_true(all(sapply(esti.metric, function(metric){
    eig.val <- eigen(metric, symmetric = T, only.values = T)
    all(eig.val$values >= 0)
  })))
  # check estimation accuracy, not so well, honestly.
  expect_equal(esti.metric, true.metric.binary, tolerance = 0.35)

  ### comparing edges
  set.seed(1)
  comp.edge <- genCompareGraph(c(500, 1000), all.edge)
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
  # arr.prob <- sigmoid.f(do.call(`-`, sq.dist))
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
      local.reach = 0.25
      , n.local = 10^5, method.trim = 'proximity'
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
  # target <- manifold$genPnt(10)
  # true.metric.compare <- apply(target, 1, manifold$metric, simplify = F)
  true.metric.compare <- lapply(true.metric, function(x) 25 * x)
  # true.metric.compare <- lapply(true.metric, function(x) x)
  esti.metric <- estiMetric(target, fit.compare)
  names(esti.metric) <- NULL
  # check symmetric and positive definite
  expect_true(all(sapply(esti.metric, isSymmetric)))
  expect_true(all(sapply(esti.metric, function(metric){
    eig.val <- eigen(metric, symmetric = T, only.values = T)
    all(eig.val$values >= 0)
  })))
  # check estimation accuracy, not so well, honestly.
  expect_equal(esti.metric, true.metric.compare, tolerance = 0.25)

  ### test additional optns
  expect_silent({
    res.more.optns <- estiMetric(target[1, ], fit.compare, optns = list(
      summary.fit = summary, tensor.only = F
    ))
  })
  expect_true(
    inherits(res.more.optns$fit[[1]], 'summary.glm')
  )

})

test_that("locMetric", {

  d <- 2

  ## Euclidean space
  set.seed(1)
  # coord of points
  coord.pnts <- matrix(runif(1000 * d, 0, 2), ncol = d)
  # pairing points
  idx.dist <- t(utils::combn(nrow(coord.pnts), 2))
  idx.dist <- idx.dist[
    sample(nrow(idx.dist), size = 50000),
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

  # # pretty much the same
  # microbenchmark::microbenchmark(list = alist(
  #   legacy = locMetric_old(arr.dist ^ 2, coord.diff = coord.diff),
  #   new = locMetric(arr.dist ^ 2, ls.cd = list(cd.cd = coord.diff))
  # ), control = list(warmup = 500L))

  # noiseless observation
  expect_equal(
    diag(rep(1, d)),
    locMetric(arr.dist ^ 2, ls.cd = list(cd.cd = coord.diff), list(
      constructor = list(metric = termMetric(d = d))
    ))$metric
  )
  # normal noise (after squared)
  expect_equal(
    diag(rep(1, d)),
    locMetric(
      arr.dist ^ 2 + rnorm(length(arr.dist), sd = sd(arr.dist) / 50),
      list(cd.cd = coord.diff), list(
        constructor = list(metric = termMetric(d = d))
    ))$metric
    , tolerance = 5e-3
  )
  # binary thresholding ###############
  expect_equal(
    diag(rep(1, d)),
    locMetric(
      arr.thre, list(cd.cd = coord.diff)
      , optns = list(
        type = 't', intercept = FALSE,
        constructor = list(metric = termMetric(d = d))
      )
    )$metric
    , tolerance = 0.07 #0.275
  )
  # binary thresholding w/ intercept ###################
  expect_equal(
    diag(rep(1, d)),
    locMetric(
      arr.thre.itcpt, list(cd.cd = coord.diff)
      , optns = list(
        type = 't', intercept = TRUE,
        constructor = list(metric = termMetric(d = d))
      )
    )$metric
    , tolerance = 0.065
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
  ls.cd.input.format <- lapply(seq(2), function(i) {
    list(cd.cd = ls.coord.diff[[i]])
  })
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
      arr.comp, ls.cd.input.format
      , optns = list(
        type = 'comp', intercept = FALSE,
        constructor = list(metric = termMetric(d = d))
      )
    )$metric
    , tolerance = 0.035 #0.175
  )
  # comparing distance w/ intercept
  res.itcpt <- locMetric(
    arr.comp, ls.cd.input.format
    , optns = list(
      type = 'comp', intercept = TRUE,
      constructor = list(metric = termMetric(d = d))
    )
  )
  expect_equal(
    diag(rep(1, d)),
    res.itcpt$metric
    , tolerance = 0.035 #0.175
  )
  # the fitted intercept should be close to 0
  expect_equal(
    0, as.numeric(res.itcpt$intercept)
    , tolerance = 0.025 #0.007
  )

  # expect_true(FALSE)
})

##### some legacy tests ########################################################

test_that("LEGACY locMetric", {

  skip_if(TRUE, message = 'Legacy tests, kept only for dev.')

  d <- 2

  ## Euclidean space
  set.seed(1)
  # coord of points
  coord.pnts <- matrix(runif(1000 * d, 0, 2), ncol = d)
  # pairing points
  idx.dist <- t(utils::combn(nrow(coord.pnts), 2))
  idx.dist <- idx.dist[
    sample(nrow(idx.dist), size = 50000),
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
    locMetric_old(arr.dist ^ 2, coord.diff)$mat.g
  )
  # normal noise (after squared)
  expect_equal(
    diag(rep(1, d)),
    locMetric_old(
      arr.dist ^ 2 + rnorm(length(arr.dist), sd = sd(arr.dist) / 50),
      coord.diff
    )$mat.g
    , tolerance = 5e-3
  )
  # binary thresholding ###############
  expect_equal(
    diag(rep(1, d)),
    locMetric_old(
      arr.thre, coord.diff
      , optns = list(type = 't', intercept = FALSE)
    )$mat.g
    , tolerance = 0.07 #0.275
  )
  # binary thresholding w/ intercept ###################
  expect_equal(
    diag(rep(1, d)),
    locMetric_old(
      arr.thre.itcpt, coord.diff
      , optns = list(type = 't', intercept = TRUE)
    )$mat.g
    , tolerance = 0.065
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
    locMetric_old(
      arr.comp, ls.coord.diff
      , optns = list(type = 'comp', intercept = FALSE)
    )$mat.g
    , tolerance = 0.035 #0.175
  )
  # comparing distance w/ intercept
  res.itcpt <- locMetric_old(
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
