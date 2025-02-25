#' @title Get weather data from SILO
#'
#' @description Extract weather data for Australia from the SILO weather data resource (https://www.longpaddock.qld.gov.au/silo/)
#' for a set of environments with defined latitude and longitude coordinates
#'
#' @param Envs Vector of environment names character strings.
#' @param Lats Vector of latitude numeric values for each environment.
#' @param Lons Vector of longitude numeric values for each environment.
#' @param Years Vector of year integer values for each environment.
#' @param vars Vector of weather variable names to get. Default is all.
#' Options include:
#'
#'    *daily_rain - Daily rainfall (mm)
#'
#'    *max_temp - Maximum temperature (deg C)
#'
#'    *min_temp - Minimum temperature (deg C)
#'
#'    *vp_deficit - Vapour pressure deficit (hPa)
#'
#'    *radiation - Solar exposure, consisting of both direct and diffuse components (MJ/m2)
#'
#' @returns A list of weather data for each weather variable and a vector of units for each variable.
#' The data list contains a a matrix or data values with environments as rows and days of the year as columns.
#'
#' @examples
#' data("CAIGE23_24envs")
#' wthr <- Get.SILO.weather()
#'
#' @references
#' Jeffrey, S. J., Carter, J. O., Moodie, K. B., & Beswick, A. R. (2001).
#'   Using spatial interpolation to construct a comprehensive archive of Australian climate data. Environmental Modelling & Software, 16(4), 309–330.
#'   https://doi.org/10.1016/S1364-8152(01)00008-1
#'
#' @export

Get.SILO.weather <- function(Envs, Lats, Lons, Years, vars = c("daily_rain", "max_temp", "min_temp", "vp_deficit", "radiation")) {
  Years <- as.integer(as.numeric(Years))
  years <- unique(Years)
  var.units <- c()
  all.vars.weather <- list()
  for (v in seq_along(vars)) {
    print(paste("Starting", vars[v]))
    all.yrs.weather <- matrix(NA, nrow = length(Envs), ncol = 365, dimnames = list(Envs, 1:365))

    for (y in seq_along(years)) {
      cat(years[y], "|", sep = "")
      addrs <- paste("https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/", vars[v], "/", years[y], ".", vars[v], ".nc", sep = "")
      nc.rast <- terra::rast(addrs)
      xvals <- terra::xFromCol(nc.rast)
      yvals <- terra::yFromRow(nc.rast)
      var.units[v] <- terra::units(nc.rast)[1]
      days <- terra::time(nc.rast)
      nc.data <- terra::as.array(nc.rast)
      dimnames(nc.data) <- list(yvals, xvals, as.character(days))
      env.info.yr.sub <- data.frame(
        "Environment" = Envs[Years == years[y]],
        "Lat" = Lats[Years == years[y]],
        "Lon" = Lons[Years == years[y]]
      )
      env.info.yr.sub$lat.ind <- sapply(env.info.yr.sub$Lat, function(x) which.min(abs(as.numeric(dimnames(nc.data)[[1]]) - as.numeric(x))))
      env.info.yr.sub$lon.ind <- sapply(env.info.yr.sub$Lon, function(x) which.min(abs(as.numeric(dimnames(nc.data)[[2]]) - as.numeric(x))))
      env.weather <- sapply(seq_len(nrow(env.info.yr.sub)), function(x) nc.data[env.info.yr.sub$lat.ind[x], env.info.yr.sub$lon.ind[x], ])
      env.weather <- t(env.weather)
      rownames(env.weather) <- env.info.yr.sub$Environment
      env.weather <- as.data.frame(env.weather)[, 1:365]
      if (ncol(env.weather) < 365) {
        env.weather <- cbind(env.weather, matrix(NA, nrow = nrow(env.weather), ncol = 365 - ncol(env.weather)))
      }
      all.yrs.weather[rownames(env.weather), ] <- as.matrix(env.weather)
    }
    gc(full = T)
    all.vars.weather[[v]] <- all.yrs.weather
    print(":)")
  }
  names(all.vars.weather) <- vars
  names(var.units) <- vars
  out <- list("data" = all.vars.weather, "units" = var.units, "Env.info" = env.info)
  return(out)
}
