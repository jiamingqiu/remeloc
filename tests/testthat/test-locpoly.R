# Functions for local polynomial

test_that("getKernelFunc", {
  arr.x <- seq(-5, 5, 0.01)
  kernel.f <- getKernelFunc("gauss", 2)
  expect_identical(
    kernel.f(arr.x),
    dnorm(arr.x, sd = 2)
  )
  expect_equal(
    integrate(kernel.f, -Inf, Inf)$value, 1
  )

  kernel.f <- getKernelFunc("step", 2)
  expect_identical(
    kernel.f(arr.x),
    0.5 * as.numeric(abs(arr.x / 2) <= 1) / 2
  )
  expect_equal(
    integrate(kernel.f, -Inf, Inf)$value, 1
  )

  # expect_true(FALSE)
})

test_that("locPolyDefine", {
  # local(browser())
  df <- data.frame(
    x1 = c(1, 2, 4),
    x2 = c(4, 5, 6),
    y = rep(0, 3)
  )
  loc.poly <- locPolyDefine(y ~ x1, df, list(
    kernel.f = getKernelFunc("step", 1)
  ))
  # matching by names
  expect_identical(
    loc.poly$mdl.mat.f(0),
    loc.poly$mdl.mat.f(c(x1 = 0, x2 = 100))
  )
  mdl.mat <- loc.poly$mdl.mat.f(0)
  mdl.mat <- matrix(as.numeric(mdl.mat), nrow = nrow(df))
  expect_identical(
    mdl.mat,
    matrix(c(
      rep(1, nrow(df)), df$x1
    ), ncol = 2)
  )
  mdl.mat <- loc.poly$mdl.mat.f(1)
  mdl.mat <- matrix(as.numeric(mdl.mat), nrow = nrow(df))
  expect_identical(
    mdl.mat,
    matrix(c(
      rep(1, nrow(df)), df$x1 - 1
    ), ncol = 2)
  )
  # the obsv weights
  expect_identical(
    loc.poly$weight.f(1.5),
    loc.poly$weight.f(c(x1 = 1.5))
  )
  expect_identical(
    loc.poly$weight.f(1.5),
    c(0.5, 0.5, 0)
  )

  # 2nd degree w/ no intercept
  loc.poly <- locPolyDefine(y ~ 0 + x1*x2, df, list(
    kernel.f = getKernelFunc("gauss", 1)
    , zero.weight.eps = 0.005
  ))
  mdl.mat <- matrix(
    as.numeric(loc.poly$mdl.mat.f(rep(0, 2))),
    nrow = nrow(df)
  )
  expect_identical(
    mdl.mat,
    with(df, matrix(c(x1, x2, x1*x2), nrow = nrow(df)))
  )
  # weights
  expect_equal(
    loc.poly$weight.f(c(x1 = 1.5, x2 = 4.5)),
    dnorm(c(sqrt(0.5), sqrt(0.5), sqrt(8.5)))
  )
  tm <- sqrt(-2 * log(sqrt(2*pi) * 0.005))
  # dnorm(tm) == 0.005
  tm <- sqrt((tm + .Machine$double.eps^(1/3)) ^ 2 / 2)
  expect_identical(
    loc.poly$weight.f(c(
      x1 = 1 + tm, x2 = 4 + tm
    ))[1]
    , 0
  )
  # stand along weighting
  side.wt.f <- locPolyWeight(
    loc.poly, getKernelFunc("gauss"), list(
      zero.weight.eps = 0.05
    )
  )
  expect_identical(
    side.wt.f(c(x1 = 1 + tm, x2 = 4 + tm)),
    loc.poly$weight.f(c(x1 = 1 + tm, x2 = 4 + tm))
  )
})
