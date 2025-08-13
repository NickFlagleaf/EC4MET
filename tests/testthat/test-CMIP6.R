test_that("get CMIP6 data works", {
  cmip.tst <- get.CMIP6.weather(Envs = c("test1", "test2"),
                                GCMs = "CNRM-ESM2-1",
                                SSPs = "ssp370",
                                Lats = c(-31.2, -31.7),
                                Lons = c(123.5, 150.5),
                                Years = c(1990, 2090),
                                ncores = 1, verbose = T,dlprompt = F)
                                  
   expect_equal(
    round(c(cmip.tst$`CNRM-ESM2-1`$ssp370$data$max_temp[1:4,10]), 2),
    c("1990_123.5_-31.2"=28.28,
      "1990_150.5_-31.7"=26.83, 
      "2090_123.5_-31.2"=45.90, 
      "2090_150.5_-31.7"=30.81)
  )
})
