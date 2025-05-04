test_that("smip function works", {
  
  obs.wthr <- get.SILO.weather(
    Envs = c("test1", "test2"),
    Lats = c(-31.2, -31.7),
    Lons = c(123.5, 150.5),
    Years = c(2020, 2021),
    ncores = 1, verbose = F
  )
  obs.wthr<-add.SMI(obs.wthr)
  
  obs.wthr$data$smi[,10]
  expect_equal(
    round(obs.wthr$data$smi[,10],2),
    c("test1"=0.14,"test2"=0.59)
  )
})


