#' @title Get weather data from BARRA-R2
#'
#' @description Extract weather data for Coupled Model Intercomparison Project phase 6 using a Quantile Delta Change method for Australia 
#' [CMIP6 QDM](https://doi.org/10.25919/03by-9y62) weather data resource for a set of environments with defined latitude and longitude coordinates. 
#' 
#' Weather variables include:
#' * `daily_rain` - Daily rainfall (mm)
#' * `max_temp` - Maximum temperature (°C)
#' * `min_temp` - Minimum temperature (°C)
#' * `vp_deficit` - Vapour pressure deficit (hPa)
#' * `radiation` - Solar exposure, consisting of both direct and diffuse components (MJ m<sup>-2</sup>)
#'
#'
#' @param Envs Vector of environment names character strings.
#' @param Lats Vector of latitude numeric values for each environment in the same order as `Envs`.
#' @param Lons Vector of longitude numeric values for each environment in the same order as `Envs`.
#' @param Years Vector of year integer values for each environment in the same order as `Envs`.
#' @param GCMs Vector of GCM names to get data from. For options see Details below.
#' @param SSPs Vector of SSP names to get data from. For options see Details below.
#' @param ncores Number (integer) of cores to use for parallel processing of gridded data up to 5 cores. Use `1` to run in series. The default (`NULL`) will
#' use the maximum available cores up to 5. If running in parallel, an output log text file will be created in the working directory.
#' @param verbose Logical.Should progress be printed?
#'
#' @details
#' The CMPI6 QDM dataset includes climate projections for two time future periods of years:
#' 1. 2035-2064
#' 2. 2070-2099
#' `Years` values must be within either of these two time periods.
#' 
#' Possible options for Global Climate Models (GCM):
#' * `ACCESS-CM2` - A much hotter future, and drier in most regions except the southeast
#' * `ACCESS-ESM1.5` - A hotter and much drier future
#' * `CESM2` - A hotter future, wetter in parts of the east and north
#' * `CMCC-ESM2` - A much hotter future with little change in mean rainfall (with regional exceptions)
#' * `CNRMESM2-1` - A much hotter future, much drier especially in the east, but wetter in the northwest
#' * `EC-Earth3` - A hotter and much wetter future for much of Australia (except southwest Western Australia)  
#' * `MPI-ESM1-2-HR` - Lower warming, mid-range changes in rainfall
#' * `NorESM2-MM` - Lower warming, mid-range changes in rainfall
#' * `UKESM1-0-LL` - Low probability, high impact case (high climate sensitivity, high Australian warming)
#' For further details see [Grose et al. 2023](https://doi.org/10.1016/j.cliser.2023.100368)
#' 
#' Possible options for Shared Socio-economic Pathways (SSP) and equivalent Representative Concentration Pathways (RCP):
#' * `ssp126` - Sustainability (RCP2.6)
#' * `ssp245` - Middle of the road (RCP4.5)
#' * `ssp370` - Regional rivalry (No equivalent RCP)
#' * `ssp585` - Fossil fuelled development (RCP8.5) 
#' 
#' For further details see [Riahi et al. 2017](https://doi.org/10.1016/j.gloenvcha.2016.05.009)
#' 
#' 
#' @returns A list of length 2:
#' * `$data` is a list of matrices of weather data for each weather variable.
#' Each data matrix has environment names as rows and days of the ywar as columns
#' * `$Env.info` is a data frame of environment names and coordinate values for environments included in the data.
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


#' @export

get.CMIP6.weather <- function(Envs,
                              Lats,
                              Lons,
                              Years,
                              GCMs,
                              SSPs,
                              ncores = NULL,
                              verbose = TRUE) {
  
  pos.GCM <- c("ACCESS-CM2","ACCESS-ESM1-5","CESM2","CMCC-ESM2","CNRM-ESM2-1","EC-Earth3","MPI-ESM1-2-HR","NorESM2-MM","UKESM1-0-LL")
  runs.codes<-c("r4i1p1f1","r6i1p1f1","r11i1p1f1","r1i1p1f1","r1i1p1f2","r1i1p1f1","r1i1p1f1","r1i1p1f1","r1i1p1f2")
  names(runs.codes)<-pos.GCM
  pos.SSP <- c("ssp585","ssp370","ssp245","ssp126")
  year.spans<-list("2035-2064" = 2035:2064,
                   "2070-2099" = 2070:2099)
  
  Years <- as.integer(as.character(Years))
  years <- unique(Years)
  Envs <- as.character(Envs)

  vars <- c("pr", "tasmax", "tasmin", "hurs", "rsds")
  method<-c("qdc-multiplicative-monthly-q100-linear-maxaf5-annual-change-matched",
              "qdc-additive-monthly-q100-linear",
              "qdc-additive-monthly-q100-linear",
              "qdc-multiplicative-monthly-q100-linear",
              "qdc-multiplicative-monthly-q100-linear-rsdscs-clipped")
  names(method)<-vars
  
  if (is.null(ncores)) {
    ncores <- min(parallel::detectCores(), length(vars))
  }
  ncores <- min(ncores, length(vars))

  if (verbose & sum(!GCMs%in%pos.GCM)>0) {
    cat("\nInvalid GCM names:",GCMs[!GCMs%in%pos.GCM])
    cat("\nPossible options include:",pos.GCM)
  }
  if (verbose & sum(!SSPs%in%pos.SSP)>0) {
    cat("\nInvalid SSP names:",SSPs[!SSPs%in%pos.SSP])
    cat("\nPossible options include:",pos.SSP)
  }
  if (verbose & !is.numeric(Lons)) {
    cat("\nLon values not numeric")
  }
  if (verbose & sum(Lons < 88.48 | Lons > 207.39) > 0) {
    cat("\nLon out of range of SILO data: 88.48 to 207.39")
  }
  if (verbose & sum(Lats < -57.97 | Lats > 12.98) > 0) {
    cat("\nLats out of range of SILO data: -57.97 to -12.98")
  }


  all.GCM.weather<-list()
    for(g in seq_along(GCMs)){
      if (verbose) {
        cat("\n###Starting ", GCMs[g],"###",sep = "")
      }
      all.SSP.weather<-list()
      for(s in seq_along(SSPs)){
        if (verbose) {
          cat("\n##Starting ", SSPs[s],"##",sep = "")
        }
        all.vars.weather <- list()
        for (v in seq_along(vars)) {
          if (verbose) {
            cat("\n#Starting ", vars[v],"#",sep = "")
          }
          
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
          
          all.yrs.weather <- foreach::foreach(y = seq_along(years), .combine = rbind, .multicombine = T) %dopar% {
            if (verbose) {
              cat(years[y], "|", sep = "")
            }
            yr.span<-names(year.spans)[which(sapply(year.spans,function(x) years[y]%in%x))]
            addrs<-paste("https://data-cbr.csiro.au/thredds/fileServer/catch_all/qdc-cmip6/QDC-CMIP6/BARRA-R2/",
                         GCMs[g],"/",SSPs[s],"/",runs.codes[GCMs[g]],"/day/",vars[v],"/AUS-05i/",yr.span,"/v20241104/",
                         vars[v],"_day_",GCMs[g],"_",SSPs[s],"_",runs.codes[GCMs[g]],"_AUS-05i_",years[y],"_",method[vars[v]],"_BARRA-R2-baseline-1985-2014_model-baseline-1985-2014.nc",sep="")
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
            env.info.yr.sub <- data.frame(
              "Environment" = Envs[Years == years[y]],
              "Lat" = Lats[Years == years[y]],
              "Lon" = Lons[Years == years[y]]
            )
            env.info.yr.sub$lon.ind <- sapply(env.info.yr.sub$Lon, function(x) which.min(abs(as.numeric(dimnames(nc.data)[[1]]) - as.numeric(x))))
            env.info.yr.sub$lat.ind <- sapply(env.info.yr.sub$Lat, function(x) which.min(abs(as.numeric(dimnames(nc.data)[[2]]) - as.numeric(x))))
            env.weather <- t(sapply(seq_len(nrow(env.info.yr.sub)), function(x) nc.data[env.info.yr.sub$lon.ind[x], env.info.yr.sub$lat.ind[x], ]))
            rownames(env.weather) <- env.info.yr.sub$Environment
            env.weather<-env.weather[,1:365]
            if (verbose & sum(is.na(env.weather)) > 0) {
              NAenvs <- Envs[!complete.cases(env.weather)]
              cat("\nNAs returned at ", paste(NAenvs, collapse = " "))
            }
            ncdf4::nc_close(ncfile)
            rm(ncfile, nc.data)
            file.remove(tmp.dir)
            return(env.weather)
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
          
        all.vars.weather[[v]] <- all.yrs.weather
        }
        names(all.vars.weather)<-vars
        all.vars.weather$rsds<-all.vars.weather$rsds/ 41.67 * 3.6 # convert from w/m2 to MJ/m2
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
        all.SSP.weather[[s]]<-out
      }
      names(all.SSP.weather)<-SSPs
      all.GCM.weather[[g]]<-all.SSP.weather
    }
  names(all.GCM.weather)<-GCMs
  gc(full = T)
  rm(all.yrs.weather,all.SSP.weather)
  
  return(all.GCM.weather)
}
