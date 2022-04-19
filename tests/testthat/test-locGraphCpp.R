test_that("Rcpp functions", {
  set.seed(1)
  x <- rnorm(5)
  expect_identical(object = orderCpp(x), expected = order(x))
  expect_identical(object = orderStable(x), expected = order(x))
  microbenchmark::microbenchmark(list = alist(
    order = order(x), orderCpp = orderCpp(x), orderStable = orderStable(x)
  ))
  x
  sort(x, index.return = T)
  match(sort(x), x)
  match(x, sort(x))
  order(x)
  rank(x)

  x[c(5, 1)]
  microbenchmark::microbenchmark(rank(x), order(order(x)))

  set.seed(1)
  coord <- matrix(runif(10), ncol = 2)
  plot(coord)
  locWindow_cpp(
    coord, rep(0.5, 2), apply(coord, 2, order), apply(coord, 2, rank)
  )
  coord <- matrix(0, nrow = 5, ncol = 2)
  coord[-1, ] <- as.matrix(expand.grid(v1 = c(-1, 1), v2 = c(-1, 1)))
  # locWindow_cpp(
  #   coord, rep(1, 2), apply(coord, 2, order), apply(coord, 2, rank)
  # )
  plot(coord)

  allEdge(coord, 1)

  set.seed(1)
  coord <- matrix(runif(2 * 10^4), ncol = 2)
  plot(coord)
  microbenchmark::microbenchmark(allEdge(coord, 0.1), times = 10L) # ~ 3s per


})
