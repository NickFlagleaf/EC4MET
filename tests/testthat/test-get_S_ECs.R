test_that("get SLGA data works", {
  testSecs <- get.S.ECs(
    Envs = c("test1", "test2", "test3"),
    Lats = c(-34.14576, -32.36873, -27.56153),
    Lons = c(138.4150, 117.8001, 152.2791),
    verbose = F, ncores = 1
  )
  expect_equal(round(c(testSecs)[c(1, 5, 20, 200)], 2), c(1.39, 1.49, 4.00, 4.30))
})
