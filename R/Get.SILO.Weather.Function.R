#' @title Get weather data from SILO
#'
#' @description Extract weather data for Australia from the [SILO](https://www.longpaddock.qld.gov.au/silo/) weather data resource
#' for a set of environments with defined latitude and longitude coordinates.
#'
#' Weather variables include:
#' * `daily_rain` - Daily rainfall (mm)
#' * `max_temp` - Maximum temperature (°C)
#' * `min_temp` - Minimum temperature (°C)
#' * `vp_deficit` - Vapour pressure deficit (hPa)
#' * `radiation` - Solar exposure, consisting of both direct and diffuse components (MJ m<sup>-2</sup>)
#'
#' @param Envs Vector of environment names character strings.
#' @param Lats Vector of latitude numeric values for each environment in the same order as `Envs`.
#' @param Lons Vector of longitude numeric values for each environment in the same order as `Envs`.
#' @param Years Vector of year integer values for each environment in the same order as `Envs`.
#' @param ncores Number (integer) of cores to use for parallel processing of gridded data up to 5 cores. Use `1` to run in series. The default (`NULL`) will
#' use the maximum available cores up to 5. If running in parallel, an output log text file will be created in the working directory.
#' @param verbose Logical. Should progress be printed?
#'
#' @details When there are only a few environments, point data will be sequentially downloaded from SILO. When there are many environments in each year,
#' data will be downloaded and extracted from whole gridded data files more efficiently. Any locations outside of the Australian land area will return `NA`.
#'
#' @returns A list of length 2:
#' * `$data` is a list of matrices of weather data for each weather variable.
#' Each data matrix has environment names as rows and days of the year as columns
#' * `$Env.info` is a data frame of environment names and coordinate values for environments included in the data.
#
#' @seealso [get.BARRA.weather()]
#'
#' @references
#' Jeffrey, S. J., Carter, J. O., Moodie, K. B., & Beswick, A. R. (2001).
#'   [Using spatial interpolation to construct a comprehensive archive of Australian climate data. Environmental Modelling & Software](https://doi.org/10.1016/S1364-8152(01)00008-1),
#'    16(4), 309–330.
#'
#' @author Nick Fradgley
#'
#' @export

get.SILO.weather <- function(Envs,
                             Lats,
                             Lons,
                             Years,
                             ncores = NULL,
                             verbose = TRUE) {
  Years <- as.integer(as.character(Years))
  years <- unique(Years)
  Envs <- as.character(Envs)
  all.vars.weather <- list()

  if (verbose & sum(!years %in% 1889:as.numeric(stringr::str_sub(Sys.Date(), 1, 4))) > 0) {
    cat("\nYears out of range of SILO data (1889 to yesterday)")
  }
  if (verbose & !is.numeric(Lats)) {
    cat("\nLat values not numeric")
  }
  if (verbose & !is.numeric(Lons)) {
    cat("\nLon values not numeric")
  }
  if (verbose & sum(Lons < 112 | Lons > 154) > 0) {
    cat("\nLon out of range of SILO data: 112 to 154")
  }
  if (verbose & sum(Lats < -44 | Lats > -10) > 0) {
    cat("\nLats out of range of SILO data: -44 to -10")
  }

  vars <- c("daily_rain", "max_temp", "min_temp", "vp_deficit", "radiation")

  if (length(Envs) * 2 < length(years) * 3 * 60) {
    if (verbose) {
      cat("\nDownloading SILO point data\n")
    }
    
    all.env.weather <- list()
    for (e in 1:length(Envs)) {
      if (verbose == TRUE & e %in% round(seq(1, length(Envs), length.out = 100))) {
        cat("|", round(e / length(Envs) * 100), "%", sep = "")
      }
      url <- paste("https://www.longpaddock.qld.gov.au/cgi-bin/silo/DataDrillDataset.php?lat=", Lats[e], "&lon=", Lons[e], "&start=", Years[e], "0101&finish=", Years[e],
        "1231&format=csv&comment=RXNDJ&username=xxx&password=apirequest",
        sep = ""
      )
      pnt.data <- utils::read.csv(url)
      all.env.weather[[e]] <- pnt.data

      if (verbose & sum(is.na(pnt.data[, vars])) > 0) {
        NAenvs <- Envs[!complete.cases(pnt.data[, vars])]
        cat("\nNAs returned at ", paste(NAenvs, collapse = " "))
      }
    }
    names(all.env.weather) <- Envs
    all.vars.weather <- lapply(vars, function(v) t(sapply(all.env.weather, function(e) e[1:365, v])))
    names(all.vars.weather) <- vars
    if (verbose) {
      cat(":)")
    }
  }

  if (length(Envs) * 2 < length(years) * 3 * 60) {
    if (verbose) {
      cat("\nDownloading SILO gridded data\n")
    }
  
    if (is.null(ncores)) {
      ncores <- min(parallel::detectCores(), length(vars))
    }
    ncores <- min(ncores, length(vars))

    if (isTRUE(ncores > 1)) { # Run in parallel
      if (verbose) {
        cat("\nRunning in parallel...")
      }
      cl <- parallel::makeCluster(ncores, outfile = "CMIP6_download_log.txt")
      doParallel::registerDoParallel(cl)
      if (verbose) {
        cat(paste("\nProgress log output to:", getwd(), "/CMIP6_download_log.txt", sep = ""))
      }
      `%dopar%` <- foreach::`%dopar%`
    }
    
    if (isTRUE(ncores == 1)) { # Run in series
      if (verbose) {
        cat("\nRunning in series")
        `%dopar%` <- foreach::`%do%`
      }
    }
      all.vars.weather <- foreach::foreach(v = seq_along(vars), .combine = list, .multicombine = T) %dopar% {
        all.yrs.weather <- matrix(NA, nrow = length(Envs), ncol = 365, dimnames = list(Envs, 1:365))
        if (verbose) {
          print(paste("Starting", vars[v]))
        }
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
          if (verbose & sum(is.na(env.weather)) > 0) {
            NAenvs <- env.info.yr.sub$Environment[!complete.cases(env.weather)]
            cat("\nNAs returned at ", paste(NAenvs, collapse = " "))
          }
          gc()
        }
        if (verbose) {
          print(":)")
        }
        return(all.yrs.weather)
      }
      
      if (isTRUE(ncores > 1)) { # if running in parallel
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
        if (verbose) {
          cat("\nFinished parallel :)")
        }
        Sys.sleep(2)
        file.remove("CMIP6_download_log.txt")
        gc(full = T)
      }
      }


  NAnums <- sapply(all.vars.weather, function(x) sum(is.na(x)))
  if (verbose & sum(NAnums) > 0) {
    cat("\nNAs in:\n")
    print(NAnums)
  }

  env.info <- data.frame("Environment" = Envs, "Lat" = Lats, "Lon" = Lons)
  out <- list("data" = all.vars.weather, "Env.info" = env.info)
  return(out)
}
