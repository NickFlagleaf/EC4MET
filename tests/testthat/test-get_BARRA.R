test_that("get BARRA data works", {
  barra.tst <- get.BARRA.weather(
    Envs = c("test1", "test2"),
    Lats = c(-31.2, -31.7),
    Lons = c(123.5, 150.5),
    Years = c(2020, 2020),
    ncores = 1, verbose = F
  )
  expect_equal(
    round(c(barra.tst$data$max_temp)[1:10], 2),
    c(38.33, 39.52, 45.24, 36.02, 28.02, 39.21, 26.29, 41.38, 30.27, 39.62)
  )
})
