test_that("get SILO data works", {
  obs.wthr <- get.SILO.weather(
    Envs = c("test1", "test2"),
    Lats = c(-31.2, -31.7),
    Lons = c(123.5, 150.5),
    Years = c(2020, 2021),
    ncores = 1, verbose = F
  )
  expect_equal(
    round(sapply(obs.wthr$data, function(x) x[,5]),2),
    matrix(c(0.0, 29.9, 11.1, 20.7, 31.1, 14.10,
             39.6,30.4, 14.8, 15.0, 29.0, 14.13),nrow = 2,byrow = T,
           dimnames = list(c("test1","test2"),
                           c("daily_rain", "max_temp", "min_temp", "vp_deficit","radiation", "day_length")))
  )
})


