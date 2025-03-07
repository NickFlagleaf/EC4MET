test_that("Get ECs works", {
  
  obs.wthr <- get.SILO.weather(Envs = c("test1","test2"),
                               Lats = c(-31.2,-31.7),
                               Lons = c(123.5,150.5),
                               Years = c(2020,2021))
  
  obs.weather.ECs <- get.W.ECs(weather = obs.wthr,
                               sow.dates = c("05/05/2020","20/06/2021"))
  
  expect_equal(round(unlist(obs.weather.ECs)[c(2,50,100,150)],2),
               c("Ndays_Sow2Emer2"=13.00,"Avtemp_Sgf2Egf2"=18.46,
                 "Ndays<0_Emer2Juv2"=10.00,"AveSR_Sgf2Egf2"=0.00))
})
