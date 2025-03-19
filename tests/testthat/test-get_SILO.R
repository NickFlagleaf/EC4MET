test_that("get SILO data works", {
  obs.wthr <- get.SILO.weather(
    Envs = c("test1", "test2"),
    Lats = c(-31.2, -31.7),
    Lons = c(123.5, 150.5),
    Years = c(2020, 2021),
    ncores = 1, verbose = F
  )
  expect_equal(
    c(obs.wthr$data$max_temp)[1:10],
    c(38.3, 23.0, 42.1, 26.0, 26.0, 26.6, 25.6, 28.0, 29.9, 30.4)
  )
})
