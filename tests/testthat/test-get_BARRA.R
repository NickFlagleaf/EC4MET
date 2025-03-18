test_that("get BARRA data works", {
  barra.tst <- get.BARRA.weather(
    Envs = c("test1", "test2"),
    Lats = c(-31.2, -31.7),
    Lons = c(123.5, 150.5),
    Years = c(2020, 2021),
    ncores = 1,verbose=T
  )
  expect_equal(
    round(c(barra.tst$data$max_temp)[1:10],2),
    c(38.33,23.32,45.24,25.37,28.02,26.49,26.29,27.27,30.27,28.41)
  )
})
