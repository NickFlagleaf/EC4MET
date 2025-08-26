
test_that("get.SILO.weather (point) is mocked and parsed correctly", {
  # this hits the little csv in tests/testthat/fixtures
  csv <- testthat::test_path("fixtures", "silo_small.csv")


  # this copies the .download_to approach but copies the local file
  fake_dl <- function(url, dest, method = "libcurl", quiet = TRUE, mode = "wb") {
    file.copy(csv, dest, overwrite = TRUE)
    dest
  }

  testthat::with_mocked_bindings(.download_to = fake_dl, {
    obs.wthr <- get.SILO.weather(
      Envs   = c("test1", "test2"),
      Lats   = c(-31.2, -31.7),
      Lons   = c(123.5, 150.5),
      Years  = c(2020, 2021),
      ncores = 1, verbose = FALSE, dlprompt = FALSE
    )

    # structure checks
    expect_true(is.list(obs.wthr))
    expect_true(all(c("data","Env.info") %in% names(obs.wthr)))
    expect_true(all(c("daily_rain","max_temp","min_temp","vp_deficit","radiation","day_length")
                    %in% names(obs.wthr$data)))
    expect_equal(rownames(obs.wthr$data$max_temp), c("test1","test2"))
    expect_equal(dim(obs.wthr$data$max_temp), c(2, 365))  # 2 envs x 365 days

# there are only two days in the local csv
expect_equal(obs.wthr$data$daily_rain[, 1:2],
             matrix(c(0, 2.4, 0, 2.4), nrow = 2, byrow = TRUE),
             ignore_attr = TRUE)

expect_equal(obs.wthr$data$max_temp[, 1:2],
             matrix(c(35.1, 29.0, 35.1, 29.0), nrow = 2, byrow = TRUE),
             ignore_attr = TRUE)

expect_equal(obs.wthr$data$min_temp[, 1:2],
             matrix(c(18.0, 17.2, 18.0, 17.2), nrow = 2, byrow = TRUE),
             ignore_attr = TRUE)

expect_equal(obs.wthr$data$vp_deficit[, 1:2],
             matrix(c(12.3, 10.1, 12.3, 10.1), nrow = 2, byrow = TRUE),
             ignore_attr = TRUE)

expect_equal(obs.wthr$data$radiation[, 1:2],
             matrix(c(20.5, 18.0, 20.5, 18.0), nrow = 2, byrow = TRUE),
             ignore_attr = TRUE)

    # later days should be NA
    expect_true(all(is.na(obs.wthr$data$max_temp[, 3:365])))

    # day_length exists and is numeric; don’t assert exact values (depends on latitude)
    expect_true(is.matrix(obs.wthr$data$day_length))
    expect_equal(dim(obs.wthr$data$day_length), c(2, 365)) 
    expect_true(all(is.finite(obs.wthr$data$day_length[, 1:10])))
  })
})
