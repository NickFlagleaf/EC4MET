#' @title Get weather data from BARRA-R2
#'
#' @description Extract weather data for Australia from the [BARRA-R2](https://opus.nci.org.au/spaces/NDP/pages/264241166/BOM+BARRA2+ob53) 
#' weather data resource for a set of environments with defined latitude and longitude coordinates. BARRA-R2 data runs from Jan 1979 to Sept 2024.
#'
#' Weather varuiables include:
#' * `daily_rain` - Daily rainfall (mm)
#' * `max_temp` - Maximum temperature (deg C)
#' * `min_temp` - Minimum temperature (deg C)
#' * `vp_deficit` - Vapour pressure deficit (hPa)
#' * `radiation` - Solar exposure, consisting of both direct and diffuse components (MJ/m2)
#'
#'
#' @param Envs Vector of environment names character strings.
#' @param Lats Vector of latitude numeric values for each environment in the same order as `Envs`.
#' @param Lons Vector of longitude numeric values for each environment in the same order as `Envs`.
#' @param Years Vector of year integer values for each environment in the same order as `Envs`.
#' @param verbose Logical.Should progress be printed?
#'
#' @returns A list of weather data for each weather variable and a vector of units for each variable.
#' The data list contains a a matrix or data values with environments as rows and days of the year as columns.
#'
#'
#' @references
#' Su, C.H., Dharssi, I., Le Marshall, J., Le, T., Rennie, S., Smith, A., Stassen, C., Steinle, P., Torrance, J., Wang, C. and Warren, R.A., 2022.
#'     BARRA2: Development of the next-generation Australian regional atmospheric reanalysis. Bureau of Meteorology.
#'     <http://www.bom.gov.au/research/publications/researchreports/BRR-067.pdf>
#'
#' @export

get.BARRA.weather <- function(Envs,
                              Lats,
                              Lons,
                              Years,
                              verbose = TRUE) {
  Years <- as.integer(as.numeric(Years))
  years <- unique(Years)
  var.units <- c()
  all.vars.weather <- list()
  mons <- stringr::str_pad(1:12, 2, pad = "0")
  vars <- c("pr", "tasmax", "tasmin", "hurs", "rsdt")
  for (v in seq_along(vars)) {
    if (verbose) {
      cat(paste("Starting", vars[v], "\n"))
    }
    all.yrs.weather <- matrix(NA, nrow = length(Envs), ncol = 365, dimnames = list(Envs, 1:365))

    for (y in seq_along(years)) {
      if (verbose) {
        cat(years[y], "|", sep = "")
      }
      all.mons.weather <- list()
      for (m in 1:length(mons)) {
        if (verbose) {
          cat(years[y], mons[m], "|", sep = "")
        }
        addrs <- paste("https://thredds.nci.org.au/thredds/fileServer/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/day/",
          vars[v], "/latest/", vars[v], "_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1_day_", years[y], mons[m], "-", years[y], mons[m], ".nc",
          sep = ""
        )
        tmp.dir <- tempfile()
        tmp.dir <- gsub("\\", "/", tmp.dir, fixed = T)
        utils::download.file(url = addrs, destfile = tmp.dir, method = "curl", quiet = T)
        ncfile <- ncdf4::nc_open(tmp.dir)
        nc.data <- ncdf4::ncvar_get(ncfile)
        start.day <- stringr::str_sub(string = ncfile$dim$time$units, start = nchar(ncfile$dim$time$units) - 9, end = nchar(ncfile$dim$time$units))
        dimnames(nc.data) <- list(
          ncfile$dim$lon$vals,
          ncfile$dim$lat$vals,
          as.character(as.Date(start.day) + ncfile$dim$time$vals)
        )
        all.mons.weather[[m]] <- nc.data
        ncdf4::nc_close(ncfile)
        rm(ncfile, nc.data)
        file.remove(tmp.dir)
        gc(full = T)
      }

      tnames <- unlist(lapply(all.mons.weather, function(x) dimnames(x)[[3]]))
      all.mons.weather <- abind::abind(all.mons.weather, along = 3)

      env.info.yr.sub <- data.frame(
        "Environment" = Envs[Years == years[y]],
        "Lat" = Lats[Years == years[y]],
        "Lon" = Lons[Years == years[y]]
      )
      env.info.yr.sub$lon.ind <- sapply(env.info.yr.sub$Lon, function(x) which.min(abs(as.numeric(dimnames(all.mons.weather)[[1]]) - as.numeric(x))))
      env.info.yr.sub$lat.ind <- sapply(env.info.yr.sub$Lat, function(x) which.min(abs(as.numeric(dimnames(all.mons.weather)[[2]]) - as.numeric(x))))
      env.weather <- sapply(seq_len(nrow(env.info.yr.sub)), function(x) all.mons.weather[env.info.yr.sub$lon.ind[x], env.info.yr.sub$lat.ind[x], ])
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


  all.vars.weather$pr <- all.vars.weather$pr * 86400 # convert "kg m-2 s-1" to "mm day-1"
  all.vars.weather$tasmax <- all.vars.weather$tasmax - 273.15 # convert from Kelvin to deg C
  all.vars.weather$tasmin <- all.vars.weather$tasmin - 273.15 # convert from Kelvin to deg C
  all.vars.weather$rsdt <- all.vars.weather$rsdt / 41.67 * 3.6 # convert from w/m2 to MJ/m2

  names(all.vars.weather) <- c("daily_rain", "max_temp", "min_temp", "relhumidity", "radiation")


  # Calculate VPD
  Tave <- (all.vars.weather$max_temp + all.vars.weather$min_temp) / 2
  es <- 0.6108 * (exp((17.27 * Tave) / (Tave + 237.3)))
  ea <- all.vars.weather$relhumidity / 100 * es
  VPD <- es - ea
  VPD <- VPD * 10 # Convert kpa to hpa

  all.vars.weather$vp_deficit <- VPD
  all.vars.weather <- all.vars.weather[!names(all.vars.weather) == "relhumidity"]
  env.info <- data.frame("Environment" = Envs, "Lat" = Lats, "Lon" = Lons)

  out <- list("data" = all.vars.weather, "Env.info" = env.info)
  return(out)
}
