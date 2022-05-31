
test_that("Riemann curvature tensors", {

  d <- 3
  manifold <- spaceHyperbolic(d = d, r = 5, model = 'ball')
  # manifold <- spaceHyperbolic(d = d, r = 5, model = 'half')
  metric <- manifold$metric
  christ.f <- getChristoffel(manifold$metric, d = d)
  cvt13.f <- get13Curvature(manifold$metric, d = d)
  cvt04.f <- get04Curvature(manifold$metric, d = d)

  # christ.f(rep(1, 2))
  # christ.f(rep(0.5, 2))
  # cvt13.f(rep(0.1, 2))
  # cvt04.f(rep(0.1, 2))

  target <- rep(0.1, d)

  tst.cvt04 <- array(NA, dim = rep(d, 4))
  for(i in seq(d)){
    for(j in seq(d)){
      for(k in seq(d)){
        for(l in seq(d)){
          tst.cvt04[i,j,k,l] <- sum(
            metric(target)[l, ] * cvt13.f(target)[i, j, k, ]
          )
        }
      }
    }
  }
  expect_equal(tst.cvt04, cvt04.f(target))
  # tst.cvt04
  # cvt04.f(target)

  tst.cvt13 <- array(NA, dim = rep(d, 4))
  grad.christ.f <- function(x) aperm(array(
    nloptr::nl.jacobian(x, christ.f), rep(d, 4)
  ), c(4, 1, 2, 3)) # [i, j, k, l] = \partial_i Christ[j, k, ^l]
  # grad.christ.f(target)

  for(i in seq(d)){
    for(j in seq(d)){
      for(k in seq(d)){
        for(l in seq(d)){

          arr.christ <- christ.f(target)
          jacob.christ <- grad.christ.f(target)
          tst.cvt13[i, j, k, l] <-
            jacob.christ[i, j, k, l] - jacob.christ[j, i, k, l] +
            sum(arr.christ[j, k, ] * arr.christ[i, , l]) -
            sum(arr.christ[i, k, ] * arr.christ[j, , l])

        }
      }
    }
  }
  expect_equal(tst.cvt13, cvt13.f(target))
  # tst.cvt13

  mat.met <- metric(target)
  inv.met <- MASS::ginv(mat.met)
  grad.met <- nloptr::nl.jacobian(target, metric)
  grad.met <- array(grad.met, dim = rep(d, 3)) # ..[i, j, k] = \partial_k g_ij
  grad.met <- reindex(grad.met, to = 'kij', from = 'ijk') #..[i,j,k] = \p_i g_jk
  tst.christ <- array(0, dim = rep(d, 3))
  for(i in seq(d)){ for(j in seq(d)){ for(k in seq(d)) {
    for(l in seq(d)) {
      tst.christ[i, j, k] = tst.christ[i, j, k] + 1/2 * inv.met[k, l] * (
        grad.met[i, j, l] + grad.met[j, i, l] - grad.met[l, i, j]
      )
    }
  }}}
  expect_equal(tst.christ, christ.f(target))
  # tst.christ - christ.f(target)
})

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

