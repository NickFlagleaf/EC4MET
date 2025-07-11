#' @title Get weather data from CMIP6 QDC
#'
#' @description Extract weather data from the Coupled Model Intercomparison Project phase 6 using a Quantile Delta Change method for Australia
#' [CMIP6 QDC](https://doi.org/10.25919/03by-9y62) weather data resource for a set of environments with defined latitude and longitude coordinates.
#'
#' @param Envs Vector of environment names character strings.
#' @param Lats Vector of latitude numeric values for each environment in the same order as `Envs`.
#' @param Lons Vector of longitude numeric values for each environment in the same order as `Envs`.
#' @param Years Vector of year integer values. Unlike [get.SILO.weather()] and [get.BARRA.weather()], `Years` should not be the same
#' length as `Envs`. Data for all locations in `Envs` will be extracted for all `Years`. Values must be within the possible time ranges of 1985-2014, 2035-2064, 
#' and 2070-2099.
#' @param GCMs Vector of GCM names to get data from. For options see Details below.
#' @param SSPs Vector of SSP names to get data from. For options see Details below.
#' @param ncores Number (integer) of cores to use for parallel processing of gridded data. Use `1` to run in series. The default (`NULL`) will
#' use the one less than the maximum available. If running in parallel, an output log text file will be created in the working directory.
#' @param verbose Logical. Should progress be printed? Default = TRUE.
#' @param dlprompt Logical. Should the user be prompted approve the total download size? Default = TRUE.
#'
#' @details
#' The CMPI6 QDC dataset is hosted on the [CSIRO data server](https://data-cbr.csiro.au/thredds/catalog/catch_all/qdc-cmip6/QDC-CMIP6/BARRA-R2/catalog.html)
#' includes climate projections for historical (1985-2014) and two future periods of years (2035-2064 and 2070-2099). All `Years` values must be within these time periods.
#' The QDC method benchmarks CMIP6 against the observational data from the [BARRA-R2](https://opus.nci.org.au/spaces/NDP/pages/264241166/BOM+BARRA2+ob53)
#' dataset so should be used in combination with the [get.BARRA.weather()] rather than the [get.SILO.weather()] function for observed environments that are
#' outside of the 1985-2014 time period of the CMIP6 QDC historical time period.
#'
#' Weather variables returned include:
#' * `daily_rain` - Daily rainfall (mm)
#' * `max_temp` - Maximum temperature (°C)
#' * `min_temp` - Minimum temperature (°C)
#' * `vp_deficit` - Vapour pressure deficit (hPa)
#' * `radiation` - Solar exposure, consisting of both direct and diffuse components (MJ m<sup>-2</sup>)
#' * `day_lengths` - Time between sunrise and sunset (h) not taken from CMIP6 QDC
#'
#' VPD in hPa is calculated as \eqn{ VPD = 10(es - ea) }, where \eqn{ es = 0.6108 \times \exp(\frac{17.27 \times T_{ave}}{T_{ave} + 237.3}) },
#' \eqn{T_{ave} } is the mean temperature in °C, \deqn{ ea = \frac{RH}{100} \times es }, and \eqn{RH} is the relative humidity (%).
#'
#' Possible options for Global Climate Models (GCM):
#' * `ACCESS-CM2` - A much hotter future, and drier in most regions except the southeast
#' * `ACCESS-ESM1-5` - A hotter and much drier future
#' * `CMCC-ESM2` - A much hotter future with little change in mean rainfall (with regional exceptions)
#' * `CNRM-ESM2-1` - A much hotter future, much drier especially in the east, but wetter in the northwest
#' * `EC-Earth3` - A hotter and much wetter future for much of Australia (except southwest Western Australia)
#' * `MPI-ESM1-2-HR` - Lower warming, mid-range changes in rainfall
#' * `NorESM2-MM` - Lower warming, mid-range changes in rainfall
#' * `UKESM1-0-LL` - Low probability, high impact case (high climate sensitivity, high Australian warming)
#'
#' For further details see [Grose et al. 2023](https://doi.org/10.1016/j.cliser.2023.100368)
#'
#' Possible options for Shared Socio-economic Pathways (SSP) and equivalent Representative Concentration Pathways (RCP)
#' with expected temperature increase range:
#' * `ssp126` - Sustainability (RCP2.6; 1.0-1.8°C)
#' * `ssp245` - Middle of the road (RCP4.5; 1.3-2.4°C)
#' * `ssp370` - Regional rivalry (No equivalent RCP; 2.8-4.6°C)
#' * `ssp585` - Fossil fueled development (RCP8.5; 3.3-5.7°C)
#'
#' For further details see [Riahi et al. 2017](https://doi.org/10.1016/j.gloenvcha.2016.05.009)
#'
#' An internet connection with high download speed is suggested for downloading gridded data for many environments.
#'
#' @returns A multi-level list of `$data` and `$Env.info` weather variables within SSPs within GCMs:
#' * `$data` is a list of matrices of weather data for each weather variable.
#' Each data matrix has environment names as rows and days of the year as columns
#' * `$Env.info` is a data frame of environment names and coordinate values for environments included in the data.
#'
#' E.g. `x$ACCESS-CM2$ssp585$daily_rain$data[1,1]` will return the rainfall on the Jan 1st at the first location in the
#' first year for the high emissions SSP585 for the ACCESS-CM2 GCM.
#'
#' @seealso [get.SILO.weather()], [get.BARRA.weather()]
#'
#' @references
#' * Grose, M. R., Narsey, S., Trancoso, R., Mackallah, C., Delage, F., Dowdy, A., Di Virgilio, G., Watterson, I., Dobrohotoff, P., Rashid, H. A., Rauniyar, S., Henley, B., Thatcher, M., Syktus, J., Abramowitz, G., Evans, J. P., Su, C.-H., & Takbash, A. (2023).
#'      [A CMIP6-based multi-model downscaling ensemble to underpin climate change services in Australia](https://doi.org/10.1016/j.cliser.2023.100368). Climate Services, 30, 100368.
#' * Irving, D., & Macadam, I. (2024). [Application-Ready Climate Projections from CMIP6 using the Quantile Delta Change method. CSIRO Climate Innovation Hub Technical Note 5](https://doi.org/10.25919/03by-9y62).
#'     Technical Note.
#' * Riahi, K., van Vuuren, D. P., Kriegler, E., Edmonds, J., O’Neill, B. C., Fujimori, S., Bauer, N., Calvin, K., Dellink, R., Fricko, O., Lutz, W., Popp, A., Cuaresma, J. C., Kc, S., Leimbach, M., Jiang, L., Kram, T., Rao, S., Emmerling, J., … Tavoni, M. (2017).
#'     [The Shared Socioeconomic Pathways and their energy, land use, and greenhouse gas emissions implications: An overview](https://doi.org/10.1016/j.gloenvcha.2016.05.009). Global Environmental Change, 42, 153–168.
#'
#' @author Nick Fradgley
#'
#' @export

get.CMIP6.weather <- function(Envs,
                              Lats,
                              Lons,
                              Years,
                              GCMs,
                              SSPs,
                              ncores = NULL,
                              verbose = TRUE,
                              dlprompt = FALSE) {
  pos.GCM <- c("ACCESS-CM2", "ACCESS-ESM1-5", "CMCC-ESM2", "CNRM-ESM2-1", "EC-Earth3", "MPI-ESM1-2-HR", "NorESM2-MM", "UKESM1-0-LL")
  runs.codes <- c("r4i1p1f1", "r6i1p1f1", "r1i1p1f1", "r1i1p1f2", "r1i1p1f1", "r1i1p1f1", "r1i1p1f1", "r1i1p1f2")
  names(runs.codes) <- pos.GCM
  pos.SSP <- c("ssp585", "ssp370", "ssp245", "ssp126")
  year.spans <- list(
    "1985-2014" = 1985:2014,
    "2035-2064" = 2035:2064,
    "2070-2099" = 2070:2099
  )

  Years <- as.integer(as.character(Years))
  Years <- Years[order(Years)]
  Envs <- as.character(Envs)

  vars <- c("pr", "tasmax", "tasmin", "hurs", "rsds")
  method <- c(
    "qdc-multiplicative-monthly-q100-linear-maxaf5-annual-change-matched",
    "qdc-additive-monthly-q100-linear",
    "qdc-additive-monthly-q100-linear",
    "qdc-multiplicative-monthly-q100-linear",
    "qdc-multiplicative-monthly-q100-linear-rsdscs-clipped"
  )
  names(method) <- vars

  dl.size <- 647893838 * length(GCMs) * length(SSPs) * length(vars) * length(Years)
  if (verbose) download_data(dlprompt, dl.size)

  if (is.null(ncores)) {
    ncores <- parallel::detectCores()-1
  }

  # Error checks
  if (verbose & !is.numeric(Lats)) stop("Lat values not numeric")
  if (verbose & !is.numeric(Lons)) stop("Lon values not numeric")
  if (verbose & sum(duplicated(Envs)) > 0) stop(paste("Duplicated Envs:", Envs[duplicated(Envs)]))

  if (verbose & length(unique(c(length(Envs), length(Lats), length(Lons)))) > 1) {
    print(sapply(list("Envs" = Envs, "Lats" = Lats, "Lons" = Lons), length))
    stop("Lengths of Envs, Lats or Lons differ")
  }

  yrs.in.spans <- sapply(Years, function(x) sum(sapply(year.spans, function(s) sum(x %in% s)) == 1))
  if (verbose & sum(yrs.in.spans == 0) > 0) {
    stop(paste("Years outside of CMIP6 range:", Years[yrs.in.spans == 0]))
  }

  if (verbose & sum(!GCMs %in% pos.GCM) > 0) {
    cat("\nPossible options include:", pos.GCM, "\n")
    stop(paste("\nInvalid GCM names:", paste(GCMs[!GCMs %in% pos.GCM], collapse = " ")))
  }
  if (verbose & sum(!SSPs %in% pos.SSP) > 0) {
    cat("\nPossible options include:", pos.SSP)
    stop(paste("Invalid SSP names:", SSPs[!SSPs %in% pos.SSP]))
  }
  if (verbose & !is.numeric(Lons)) {
    stop("Lon values not numeric")
  }
  if (verbose & sum(Lons < 88.48 | Lons > 207.39) > 0) {
    stop("Lons out of range of CMIP6 QDC data: 88.48 to 207.39")
  }
  if (verbose & sum(Lats < -57.97 | Lats > 12.98) > 0) {
    stop("Lats out of range of CMIP6 QD data: -57.97 to -12.98")
  }

  if(!capabilities("libcurl")){warning("libcurl is not supported!")}

  all.GCM.weather <- list()
  for (g in seq_along(GCMs)) {
    if (verbose) {
      cat("\n###Starting ", GCMs[g], "###", sep = "")
    }
    all.SSP.weather <- list()
    for (s in seq_along(SSPs)) {
      if (verbose) {
        cat("\n##Starting ", SSPs[s], "##", sep = "")
      }
      all.vars.weather <- list()
      for (v in seq_along(vars)) {
        if (verbose) {
          cat("\n#Starting ", vars[v], "#", sep = "")
        }

        batch.years <- lapply(year.spans, function(x) Years[Years %in% x])

        # Bulk download for 2035-2064
        if (length(batch.years$`1985-2014`) > 0) {
          addrs1985.2014 <- paste("https://data-cbr.csiro.au/thredds/fileServer/catch_all/qdc-cmip6/QDC-CMIP6/BARRA-R2/obs/historical/v1/day/", vars[v],
            "/AUS-05i/1985-2014/v20241104/", vars[v], "_day_BARRA-R2_historical_v1_AUS-05i_", batch.years$`1985-2014`, ".nc",
            sep = ""
          )
        } else {
          addrs1985.2014 <- NULL
        }

        if (length(batch.years$`2035-2064`) > 0) {
          addrs2035.2064 <- paste("https://data-cbr.csiro.au/thredds/fileServer/catch_all/qdc-cmip6/QDC-CMIP6/BARRA-R2/",
            GCMs[g], "/", SSPs[s], "/", runs.codes[GCMs[g]], "/day/", vars[v], "/AUS-05i/2035-2064/v20241104/",
            vars[v], "_day_", GCMs[g], "_", SSPs[s], "_", runs.codes[GCMs[g]], "_AUS-05i_", batch.years$`2035-2064`, "_", method[vars[v]],
            "_BARRA-R2-baseline-1985-2014_model-baseline-1985-2014.nc",
            sep = ""
          )
        } else {
          addrs2035.2064 <- NULL
        }

        if (length(batch.years$`2070-2099`) > 0) {
          addrs2070.2099 <- paste("https://data-cbr.csiro.au/thredds/fileServer/catch_all/qdc-cmip6/QDC-CMIP6/BARRA-R2/",
            GCMs[g], "/", SSPs[s], "/", runs.codes[GCMs[g]], "/day/", vars[v], "/AUS-05i/2070-2099/v20241104/",
            vars[v], "_day_", GCMs[g], "_", SSPs[s], "_", runs.codes[GCMs[g]], "_AUS-05i_", batch.years$`2070-2099`, "_", method[vars[v]],
            "_BARRA-R2-baseline-1985-2014_model-baseline-1985-2014.nc",
            sep = ""
          )
        } else {
          addrs2070.2099 <- NULL
        }

        addrs <- c(addrs1985.2014, addrs2035.2064, addrs2070.2099)

        tmp.dir <- tempfile()
        tmp.dir <- gsub("\\", "/", tmp.dir, fixed = T)
        tmp.dir <- paste(tmp.dir, "_", Years, sep = "")
        cat("\nDownloading .nc files...")
        options(timeout = max(80000, getOption("timeout")))

        utils::download.file(url = addrs, destfile = tmp.dir, method = "libcurl", quiet = T, mode = "wb")
        finfo<-file.info(tmp.dir)
        tryagain<-which(finfo$size<500000000 | is.na(finfo$size))
        if(length(tryagain)>0){
          utils::download.file(url = addrs[tryagain], destfile = tmp.dir[tryagain], method = "libcurl", quiet = T, mode = "wb")
        }
        finfo<-file.info(tmp.dir)
        tryagain<-which(finfo$size<500000000 | is.na(finfo$size))
        if(length(tryagain)>0){
          utils::download.file(url = addrs[tryagain], destfile = tmp.dir[tryagain], method = "libcurl", quiet = T, mode = "wb")
        }
 
        if (isTRUE(ncores > 1)) { # Run in parallel
          if (verbose) {
            cat("\nRunning in parallel...")
          }
          suppressWarnings(file.remove("CMIP6_download_log.txt"))
          cl <- parallel::makeCluster(ncores, outfile = "CMIP6_download_log.txt")
          doParallel::registerDoParallel(cl)
          if (verbose) {
            cat(paste("\nProgress log output to:", getwd(), "/CMIP6_download_log.txt", sep = ""))
          }
          on.exit(closeAllConnections())
          `%dopar%` <- foreach::`%dopar%`
        }

        if (isTRUE(ncores == 1)) { # Run in series
          if (verbose) {
            cat("\nRunning in series...")
            `%dopar%` <- foreach::`%do%`
          }
        }

        all.yrs.weather <- foreach::foreach(y = seq_along(tmp.dir), .combine = rbind, .multicombine = T, .export = "nc.process") %dopar% {
          if (verbose) {
            cat(Years[y], "|", sep = "")
          }
          nc.data <- try(nc.process(tmp.dir[y]))
          
          if(class(nc.data)=="try-error"){
            utils::download.file(url = addrs[y], destfile = tmp.dir[y], method = "libcurl", quiet = T, mode = "wb")
            nc.data <- try(nc.process(tmp.dir[y]))
          }
          
          yr <- stringr::str_sub(dimnames(nc.data)[[3]][1], 1, 4)
          file.remove(tmp.dir[y])
          lon.ind <- sapply(Lons, function(x) which.min(abs(as.numeric(dimnames(nc.data)[[1]]) - as.numeric(x))))
          lat.ind <- sapply(Lats, function(x) which.min(abs(as.numeric(dimnames(nc.data)[[2]]) - as.numeric(x))))
          env.weather <- t(sapply(seq_len(length(lon.ind)), function(x) nc.data[lon.ind[x], lat.ind[x], ]))
          rownames(env.weather) <- paste(yr, Lons, Lats, sep = "_")
          env.weather <- env.weather[, 1:365]
          if (verbose & sum(is.na(env.weather)) > 0) {
            NAenvs <- Envs[!complete.cases(env.weather)]
            cat("\nNAs returned at ", paste(NAenvs, collapse = " "))
          }
          return(env.weather)
        }

        if (isTRUE(ncores > 1)) { # if running in parallel
          parallel::stopCluster(cl)
          doParallel::stopImplicitCluster()
          if (verbose) {
            cat("\nFinished parallel :)")
          }
          Sys.sleep(5)
          file.remove("CMIP6_download_log.txt")
        }
        gc(full = T)
        colnames(all.yrs.weather) <- 1:ncol(all.yrs.weather)
        all.vars.weather[[v]] <- all.yrs.weather
      }
      names(all.vars.weather) <- vars
      all.vars.weather$rsds <- all.vars.weather$rsds / 41.67 * 3.6 # convert from w/m2 to MJ/m2
      names(all.vars.weather) <- c("daily_rain", "max_temp", "min_temp", "relhumidity", "radiation")

      # Calculate VPD
      all.vars.weather$vp_deficit <- vpdfun(
        tmin = all.vars.weather$min_temp,
        tmax = all.vars.weather$max_temp,
        relh = all.vars.weather$relhumidity
      )
      all.vars.weather <- all.vars.weather[!names(all.vars.weather) == "relhumidity"]

      Lats.full <- as.numeric(sapply(rownames(all.vars.weather$radiation), function(x) stringr::str_split(x, pattern = "_")[[1]][3]))
      DLs <- t(sapply(Lats.full, function(x) chillR::daylength(latitude = x, JDay = 1:370, notimes.as.na = FALSE)$Daylength))
      rownames(DLs) <- rownames(all.vars.weather$radiation)
      all.vars.weather$day_length <- DLs

      env.info <- data.frame(
        "Environment" = rownames(all.vars.weather$radiation),
        "Year" = as.numeric(sapply(rownames(all.vars.weather$radiation), function(x) stringr::str_split(x, pattern = "_")[[1]][1])),
        "Lon" = as.numeric(sapply(rownames(all.vars.weather$radiation), function(x) stringr::str_split(x, pattern = "_")[[1]][2])),
        "Lat" = as.numeric(sapply(rownames(all.vars.weather$radiation), function(x) stringr::str_split(x, pattern = "_")[[1]][3]))
      )
      out <- list("data" = all.vars.weather, "Env.info" = env.info)
      all.SSP.weather[[s]] <- out
    }
    names(all.SSP.weather) <- SSPs
    all.GCM.weather[[g]] <- all.SSP.weather
  }
  names(all.GCM.weather) <- GCMs
  gc(full = T)
  rm(all.yrs.weather, all.SSP.weather)

  return(all.GCM.weather)
}
