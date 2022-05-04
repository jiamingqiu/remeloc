test_that('geodesic', {

  d <- 2
  manifold <- spaceEuclidean(d)
  geo.f <- getGeodesic(manifold$metric, d = d)
  geo.curve <- geo.f(rep(0, d), rep(1, d))
  expect_identical(geo.curve[, 'x1'], geo.curve[, 'time'])
  expect_identical(geo.curve[, 'x2'], geo.curve[, 'time'])
  expect_identical(geo.curve[, 'v1'], rep(1, nrow(geo.curve)))
  expect_identical(geo.curve[, 'v2'], rep(1, nrow(geo.curve)))

})

test_that("vectorizing symmetric matrix", {

  mat <- matrix(runif(3 * 3), 3, 3)
  mat <- (mat + t(mat)) / 2
  vec <- c(diag(mat), mat[lower.tri(mat)])

  expect_identical(vec2SymMat(vec), mat)
  expect_identical(symMat2Vec(mat), vec)
  expect_identical(mat[vecQuadIdx(3, T)], vec)

})

