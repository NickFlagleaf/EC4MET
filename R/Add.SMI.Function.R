#' @title Add estimated Soil Moisture Index (SMI) to weather data
#'
#' @description This function uses daily rainfall and temperature to estimate a SMI from a weather data object that has been output from
#' get weather functions such as [get.SILO.weather()].
#'
#' @param weather A list of length 2. Weather data as outputted from the [get.SILO.weather()].
#'
#' `weather$data` is a list of data matrices for each covariate with rows as environments and days of the year as columns. Weather covariate names for list items should include:
#' * `daily_rain`
#' * `max_temp`
#'
#' `weather$Env.info` is a data frame of info for each environment but is not required for this function.
#'
#' @details
#' A simplified SMI is calculated based on a dynamic water balance model where soil moisture is proportionally dried by daily evapotranspiration 
#' and increased by daily precipitation. The equation for SMI for a following day(\eqn{SMI_{i+1}}) is:
#' \deqn{SMI_{i+1} = \min{(1,\frac{SMI_i}{\frac{PET_i}{AWC \times c_s}c_d + 1} + (P_i \times c_p))}}
#' where \eqn{SMI_i} is the SMI for day \eqn{i}, \eqn{PET_i} is the potential evapotranspiration using the Thornthwaite equation for day \eqn{i}, and
#' \eqn{P_i} is the precipitation (mm) for day \eqn{i}. The constants \eqn{c_s}, \eqn{c_d}, and \eqn{c_p} were optimised based on ground truthed data from
#' [TERN](https://www.tern.org.au/news-smips-soil-moisture/)
#'
#' @returns The same `weather` list object as inputted with the addition of a matrix of daily SMI values in `...$data$smi`.
#
#' @seealso [get.BARRA.weather()], [get.SILO.weather()], [get.CMIP6.weather()]
#'
#' @author Nick Fradgley
#'
#' @export

add.SMI <- function(weather) {
  
  awc.urls<-c(
    "/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/AWC/AWC_000_005_EV_N_P_AU_TRN_N_20210614.tif",
    "/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/AWC/AWC_005_015_EV_N_P_AU_TRN_N_20210614.tif",
    "/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/AWC/AWC_015_030_EV_N_P_AU_TRN_N_20210614.tif",
    "/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/AWC/AWC_030_060_EV_N_P_AU_TRN_N_20210614.tif",
    "/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/AWC/AWC_060_100_EV_N_P_AU_TRN_N_20210614.tif",
    "/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/AWC/AWC_100_200_EV_N_P_AU_TRN_N_20210614.tif"
  )
  
  rs <- sapply(awc.urls,function(x) terra::rast(x))
  AWCdata <- matrix(NA, nrow = nrow(weather$Env.info), ncol = 6, dimnames = list(
    weather$Env.info$Environment,
    c("0-0.05m", "0.05-0.15m", "0.15-0.3m", "0.3-0.6m", "0.6-1m", "1-2m")
  ))
  for (d in 1:length(rs)) {
    cat("\rGetting available water data: ", round(d / length(rs) * 100, 0), "%", sep = "")
    vals <- terra::extract(x = rs[[d]], y = weather$Env.info[, c("Lon", "Lat")], search_radius = 5000)
    AWCdata[, d] <- vals[, 2]
  }

  Envs <- weather$Env.info$Environment
  ndays <- 100
  dayrange <- 1:365
  smips <- t(sapply(Envs, function(e) {
    cat("\r", rep("", 50))
    cat("\rProgress: ", round(which(Envs == e) / length(Envs) * 100), "%", sep = "")
    pr <- weather$data$daily_rain[e, ]
    tmax <- weather$data$max_temp[e, ]
    tmin <- weather$data$min_temp[e, ]
    dl <- weather$data$day_length[e, ]
    awc <- mean(AWCdata[e, ])
    smi <- SMIfun(prs = pr, tmins = tmin, tmaxs = tmax, dls = dl, AWC = awc)
    return(smi)
  }))
  colnames(smips) <- dayrange
  weather$data$smi <- smips
  return(weather)
}

