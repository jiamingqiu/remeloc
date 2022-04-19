test_that("locWindow_cpp and allEdge", {

  # error handling safe guide (especially marginOrder and marginRank)
  set.seed(1)
  coord <- matrix(runif(2 * 10^3), ncol = 2)
  marginOrder = apply(coord, 2, order)
  marginRank = apply(coord, 2, rank)
  loc.reach <- rep(1, 2)
  expect_error(
    locWindow_cpp(coord, 1, marginOrder, marginRank)
  )
  expect_error(
    locWindow_cpp(coord, loc.reach, marginOrder - 1, marginRank)
  )
  expect_error(
    locWindow_cpp(coord, loc.reach, marginOrder, head(marginRank))
  )
  # if not integer, will become integer, though extremely risky
  # locWindow_cpp(
  #   coord, loc.reach, pmin(marginOrder + 0.5, nrow(coord)), marginRank
  # )

  ## checking correctness of results
  coord <- matrix(0, nrow = 5, ncol = 2)
  coord[-1, ] <- as.matrix(expand.grid(v1 = c(-1, 1), v2 = c(-1, 1)))
  # plot(coord)

  # check if ties.method = 'first' is as desired (for stable sort only)
  expect_identical(
    apply(apply(coord, 2, order), 2, order),
    apply(coord, 2, rank, ties.method = 'first')
  )

  # locWindow_cpp(
  #   coord, rep(1, 2), apply(coord, 2, order),
  #   apply(coord, 2, rank, ties.method = 'first')
  # ) %>% .[[1]] %>% class
  # the following checks results as well as result type (integer)
  expect_identical(
    allEdge(coord, 1),
    matrix(c(
      rep(1L, nrow(coord)), seq(2, 5), seq(1, 5), seq(2, 5)
    ), ncol = 2)
  )
  expect_identical(
    allEdge(coord, 0.9),
    matrix(seq(nrow(coord)), ncol = 2, nrow = nrow(coord))
  )
})

# some dev codes
# set.seed(1)
# x <- rnorm(5)
# expect_identical(object = orderCpp(x), expected = order(x))
# expect_identical(object = orderStable(x), expected = order(x))
# microbenchmark::microbenchmark(list = alist(
#   order = order(x), orderCpp = orderCpp(x), orderStable = orderStable(x)
# ))
# x
# sort(x, index.return = T)
# match(sort(x), x)
# match(x, sort(x))
# order(x)
# rank(x)
#
# x[c(5, 1)]
# microbenchmark::microbenchmark(rank(x), order(order(x)))
#
# set.seed(1)
# coord <- matrix(runif(10), ncol = 2)
# plot(coord)
# locWindow_cpp(
#   coord, rep(0.5, 2), apply(coord, 2, order), apply(coord, 2, rank)
# )
# coord <- matrix(0, nrow = 5, ncol = 2)
# coord[-1, ] <- as.matrix(expand.grid(v1 = c(-1, 1), v2 = c(-1, 1)))
# # locWindow_cpp(
# #   coord, rep(1, 2), apply(coord, 2, order), apply(coord, 2, rank)
# # )
# plot(coord)
#
# allEdge(coord, 1)
#
# set.seed(1)
# coord <- matrix(runif(2 * 10^4), ncol = 2)
# plot(coord)
# microbenchmark::microbenchmark(allEdge(coord, 0.1), times = 10L) # ~ 3s per
