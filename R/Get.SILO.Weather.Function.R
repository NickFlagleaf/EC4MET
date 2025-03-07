#' @title Get weather data from SILO
#'
#' @description Extract weather data for Australia from the [SILO](https://www.longpaddock.qld.gov.au/silo/) weather data resource
#' for a set of environments with defined latitude and longitude coordinates.
#'
#' Weather variables include:
#' * `daily_rain` - Daily rainfall (mm)
#' * `max_temp` - Maximum temperature (deg C)
#' * `min_temp` - Minimum temperature (deg C)
#' * `vp_deficit` - Vapour pressure deficit (hPa)
#' * `radiation` - Solar exposure, consisting of both direct and diffuse components (MJ/m2)
#'
#' @param Envs Vector of environment names character strings.
#' @param Lats Vector of latitude numeric values for each environment in the same order as `Envs`.
#' @param Lons Vector of longitude numeric values for each environment in the same order as `Envs`.
#' @param Years Vector of year integer values for each environment in the same order as `Envs`.
#' @param verbose Logical. Should progress be printed?
#'
#' @details When there are only a few environments, point data will be sequentially downloaded from SILO. When there are many environments in each year, 
#' data will be downloaded and extracted from whole gridded data files more efficiently.
#' 
#' 
#' 
#' 
#'
#' @returns A list of weather data for each weather variable and a vector of units for each variable.
#' The data list contains a a matrix or data values with environments as rows and days of the year as columns.
#'
#' @references
#' Jeffrey, S. J., Carter, J. O., Moodie, K. B., & Beswick, A. R. (2001).
#'   Using spatial interpolation to construct a comprehensive archive of Australian climate data. Environmental Modelling & Software, 16(4), 309–330.
#'   <https://doi.org/10.1016/S1364-8152(01)00008-1>
#'
#' @export

get.SILO.weather <- function(Envs,
                             Lats,
                             Lons,
                             Years,
                             verbose = TRUE) {
  Years <- as.integer(as.numeric(Years))
  years <- unique(Years)
  all.vars.weather <- list()

  vars <- c("daily_rain", "max_temp", "min_temp", "vp_deficit", "radiation")

  if (length(Envs)*2 < length(years)*3*60) {
    if (verbose) {
      print("Downloading SILO point data")
    }
    all.env.weather <- list()
    for (e in 1:length(Envs)) {
      if (verbose) {
        cat("|")
      }
      url <- paste("https://www.longpaddock.qld.gov.au/cgi-bin/silo/DataDrillDataset.php?lat=", Lats[e], "&lon=", Lons[e], "&start=", Years[e], "0101&finish=", Years[e],
        "1231&format=csv&comment=RXNDJ&username=xxx&password=apirequest",
        sep = ""
      )
      pnt.data <- utils::read.csv(url)
      all.env.weather[[e]] <- pnt.data
    }
    names(all.env.weather) <- Envs
    all.vars.weather <- lapply(vars, function(v) t(sapply(all.env.weather, function(e) e[1:365, v])))
    names(all.vars.weather) <- vars
    if (verbose) {
      print(":)")
    }
  }


  if (length(Envs)*2 > length(years)*3*60) {
    if (verbose) {
      print("Downloading SILO gridded data")
    }

    for (v in seq_along(vars)) {
      if (verbose) {
        print(paste("Starting", vars[v]))
      }
      all.yrs.weather <- matrix(NA, nrow = length(Envs), ncol = 365, dimnames = list(Envs, 1:365))

      for (y in seq_along(years)) {
        if (verbose) {
          cat(years[y], "|", sep = "")
        }
        addrs <- paste("https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/", vars[v], "/", years[y], ".", vars[v], ".nc", sep = "")
        nc.rast <- terra::rast(addrs)
        xvals <- terra::xFromCol(nc.rast)
        yvals <- terra::yFromRow(nc.rast)
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
      if (verbose) {
        print(":)")
      }
    }
    names(all.vars.weather) <- vars
  }

  env.info <- data.frame("Environment" = Envs, "Lat" = Lats, "Lon" = Lons)

  out <- list("data" = all.vars.weather, "Env.info" = env.info)
  return(out)
}
