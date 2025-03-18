test_that("Get W ECs works", {
  tst.wthr <- get.SILO.weather(
    Envs = c("test1", "test2"),
    Lats = c(-31.2, -31.7),
    Lons = c(123.5, 150.5),
    Years = c(2020, 2021),
    ncores = 1,verbose=F
  )

  obs.weather.ECs <- get.W.ECs(
    weather = tst.wthr,
    sow.dates = c("05/05/2020", "20/06/2021")
  )

  expect_equal(
    round(unlist(obs.weather.ECs$ECs)[c(2, 50, 10, 150)], 2),
    c(
      "Ndays_Sow2Emer2" = 13.00, "Ndd_Sgf2Egf2" = 8,
      "Ndays_Flw2Sgf2" = 13, "AveSR_Sgf2Egf2" = 18.96
    )
  )
})
