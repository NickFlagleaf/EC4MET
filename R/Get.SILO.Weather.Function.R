#' @title Get weather data from SILO
#'
#' @description Extract weather data for Australia from the [SILO](https://www.longpaddock.qld.gov.au/silo/) weather data resource
#' for a set of environments with defined latitude and longitude coordinates.
#'
#' @param Envs Vector of environment names character strings.
#' @param Lats Vector of latitude numeric values for each environment in the same order as `Envs`.
#' @param Lons Vector of longitude numeric values for each environment in the same order as `Envs`.
#' @param Years Vector of year integer values for each environment in the same order as `Envs`.
#' @param plus.yr Logical. Should the subsequent years weather data also be downloaded? This may be needed if the estimated crop growth stages in the `get.W.ECs()`
#' function extend after the same year of sowing. 
#' @param ncores Number (integer) of cores to use for parallel processing of gridded data up to 5 cores. Use `1` to run in series. The default (`NULL`) will
#' use the maximum available cores up to 5. If running in parallel, an output log text file will be created in the working directory.
#' @param verbose Logical. Should progress be printed? Default = TRUE.
#' @param dlprompt Logical. Should the user be prompted approve the total download size? Default = FALSE
#'
#' @details
#' Weather variables returned include:
#' * `daily_rain` - Daily rainfall (mm)
#' * `max_temp` - Maximum temperature (°C)
#' * `min_temp` - Minimum temperature (°C)
#' * `vp_deficit` - Vapour pressure deficit (hPa)
#' * `radiation` - Solar exposure, consisting of both direct and diffuse components (MJ m<sup>-2</sup>)
#' * `day_lengths` - Time between sunrise and sunset (h) not taken from SILO
#'
#' When there are only a few environments, point data will be sequentially downloaded from SILO. When there are many environments in each year,
#' data will be downloaded and extracted from whole gridded data files more efficiently. Any locations outside of the Australian land area will return `NA`.
#'
#' An internet connection with high download speed is suggested for downloading gridded data for many environments.
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
#'   [Using spatial interpolation to construct a comprehensive archive of Australian climate data](https://doi.org/10.1016/S1364-8152(01)00008-1).
#'   \emph{Environmental Modelling & Software}, 16(4), 309–330.
#'
#' @author Nick Fradgley
#'
#' @export

get.SILO.weather <- function(Envs,
                             Lats,
                             Lons,
                             Years,
                             plus.yr = FALSE,
                             ncores = NULL,
                             verbose = TRUE,
                             dlprompt = FALSE) {
  Years <- as.integer(as.character(Years))
  years <- unique(Years)
  Envs <- as.character(Envs)
  all.vars.weather <- list()
  vars <- c("daily_rain", "max_temp", "min_temp", "vp_deficit", "radiation")

  if (length(unique(c(length(Envs), length(Lats), length(Lons), length(Years)))) > 1) {
    print(sapply(list("Envs" = Envs, "Lats" = Lats, "Lons" = Lons, "Years" = Years), length))
    stop("Lengths of Envs, Lats, Lons or Years differ")
  }
  if (sum(!years %in% 1889:as.numeric(stringr::str_sub(Sys.Date(), 1, 4))) > 0) stop("Years out of range of SILO data (1889 to yesterday)")
  if (!is.numeric(Lats)) stop("Lat values not numeric")
  if (!is.numeric(Lons)) stop("Lon values not numeric")
  if (sum(duplicated(Envs)) > 0) stop(paste("Duplicated Envs:", Envs[duplicated(Envs)]))
  if (sum(Lons < 112 | Lons > 154) > 0) stop("Lon out of range of SILO data: 112 to 154")
  if (sum(Lats < -44 | Lats > -10) > 0) stop("Lats out of range of SILO data: -44 to -10")

  this.year<-as.numeric(stringr::str_sub(Sys.Date(),1,4))
  if(plus.yr) this.year <- this.year - 1
  if(sum(Years >= this.year)>0) cat(crayon::yellow("Downloading data from this year. May be incomplete.\n"))
  
  if (length(Envs) < 500) {
    dl.size <- 29982 * length(Envs)

    if (verbose) download_data(dlprompt, dl.size)

    if (verbose) {
      cat("\nDownloading SILO point data\n")
    }
    urls <- paste("https://www.longpaddock.qld.gov.au/cgi-bin/silo/DataDrillDataset.php?lat=", Lats, "&lon=", Lons, "&start=", Years, "0101&finish=", Years,
      "1231&format=csv&comment=RXNDJ&username=xxx&password=apirequest",
      sep = ""
    )
    
    if(plus.yr){
    yr2urls <- paste("https://www.longpaddock.qld.gov.au/cgi-bin/silo/DataDrillDataset.php?lat=", Lats, "&lon=", Lons, "&start=", Years+1, "0101&finish=", Years+1,
                            "1231&format=csv&comment=RXNDJ&username=xxx&password=apirequest",
                            sep = "")
    }
  
  
    tmp.files <- paste0(gsub("\\", "/",tempdir(), fixed = T),"/",Envs,".csv")
    if(plus.yr) tmp.files.plus <- paste0(gsub("\\", "/",tempdir(), fixed = T),"/",Envs,"_plus.csv")
    
    out<-try(.download_to(urls, tmp.files,quiet = T))
    if(plus.yr) out<-try(.download_to(yr2urls, tmp.files.plus,quiet = T))

    finfo<-file.info(tmp.files)
    tryagain<-which(finfo$size<29000  | is.na(finfo$size))
    if (length(tryagain) > 0) try(.download_to(urls[tryagain], tmp.files[tryagain]))
    
    if(plus.yr) {
    finfo2<-file.info(tmp.files.plus)
    tryagain2<-which(finfo2$size<29000  | is.na(finfo2$size))
    if (length(tryagain2) > 0) try(.download_to(yr2urls[tryagain2], tmp.files.plus[tryagain2]))
    }
    
    
    all.env.weather <- list()
    for (e in 1:length(Envs)) {
      if (verbose == TRUE & e %in% round(seq(1, length(Envs), length.out = 100))) {
        cat("\r|", round(e / length(Envs) * 100), "%", sep = "")
      }
      pnt.data <- utils::read.csv(tmp.files[e])
      if(plus.yr) pnt.data <- rbind(pnt.data,utils::read.csv(tmp.files.plus[e]))
      
      all.env.weather[[e]] <- pnt.data

      if (verbose & sum(is.na(pnt.data[, vars])) > 0) {
        NAenvs <- Envs[!stats::complete.cases(pnt.data[, vars])]
        cat("\nNAs returned at ", paste(NAenvs, collapse = " "))
      }
      rm(pnt.data)
    }
    file.remove(tmp.files)
    gc()
    names(all.env.weather) <- Envs
    if(plus.yr){ndays<-730}else{ndays<-365}
    all.vars.weather <- lapply(vars, function(v) t(sapply(all.env.weather, function(e) e[1:ndays, v])))
    names(all.vars.weather) <- vars
    if (verbose) cat(crayon::green(" :)"))
    }

  if (length(Envs) > 499) {
    dl.size <- 419290699 * length(years) * length(vars)
    if (verbose) download_data(dlprompt, dl.size)

    if (verbose) cat("\nDownloading SILO gridded data")

    if (is.null(ncores)) {
      ncores <- min(parallel::detectCores(), length(vars))
    }
    ncores <- min(ncores, length(vars))

    if (isTRUE(ncores > 1)) { # Run in parallel
      if (verbose) cat("\nRunning in parallel...")
      if ("SILO_download_log.txt" %in% dir()) suppressWarnings(file.remove("SILO_download_log.txt"))
      cl <- parallel::makeCluster(ncores, outfile = "SILO_download_log.txt")
      doParallel::registerDoParallel(cl)
      if (verbose) {
        cat(paste("\nProgress log output to: ", getwd(), "/SILO_download_log.txt", sep = ""))
      }
      `%how%` <- foreach::`%dopar%`
      parallel::clusterExport(cl, c("nc.process", ".download_to"), envir = environment())
      on.exit(closeAllConnections())
    }

    if (isTRUE(ncores == 1)) { # Run in series
      if (verbose) cat("\nRunning in series\n")
      `%how%` <- foreach::`%do%`
    }

    all.vars.weather <- foreach::foreach(
      v = seq_along(vars),
      .combine = list,
      .multicombine = T, 
      .export = c("nc.process", ".download_to")
      ) %how% {
        
        if(plus.yr){ndays<-730}else{ndays<-365}
        all.yrs.weather <- matrix(NA, nrow = length(Envs), ncol = ndays, dimnames = list(Envs, 1:ndays))
        if (verbose) cat("Starting", vars[v])
        if (verbose) cat("\nDownloading .nc files...\n")
        addrs <- paste("https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/", vars[v], "/", years, ".", vars[v], ".nc", sep = "")
        if(plus.yr) {
          addrs.plus <- paste("https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/", vars[v], "/", years+1, ".", vars[v], ".nc", sep = "")
          addrs<-unique(c(addrs,addrs.plus))
        }
        addrs<-addrs[order(addrs)]
        
        tmp.dir <- tempdir()
        tmp.dir <- gsub("\\", "/", tmp.dir, fixed = T)
        tmp.files <- paste(tmp.dir, "/SILO_", vars[v], "_", years, ".nc", sep = "")
        if(plus.yr) {
          tmp.files.plus <- paste(tmp.dir, "/SILO_", vars[v], "_", years+1, ".nc", sep = "")
          tmp.files<-unique(c(tmp.files,tmp.files.plus))
        }
        tmp.files<-tmp.files[order(tmp.files)]
        
        options(timeout = max(80000, getOption("timeout")))
        .download_to(addrs, tmp.files)
  
        finfo<-file.info(tmp.files)
        tryagain<-which(finfo$size<40000000 | is.na(finfo$size))
        if (length(tryagain) > 0) try(.download_to(addrs[tryagain], tmp.files[tryagain]))
        
        for (y in seq_along(years)) {
          if (verbose) cat(years[y], "|", sep = "")
          env.info.yr.sub <- data.frame(
            "Environment" = Envs[Years == years[y]],
            "Lat" = Lats[Years == years[y]],
            "Lon" = Lons[Years == years[y]]
          )
          nc.data <- try(nc.process(tmp.files[y]))
          if (inherits(nc.data, "try-error")) {
  +          try(.download_to(addrs[y], tmp.files[y]))
            nc.data <- try(nc.process(tmp.files[y]))
          }
          if(plus.yr){
            nc.data.plus <- try(nc.process(tmp.files[y+1]))
            if (inherits(nc.data.plus, "try-error")) {
              +          try(.download_to(addrs[y+1], tmp.files[y+1]))
              nc.data.plus <- try(nc.process(tmp.files[y+1]))
            }
            nc.data<-abind::abind(nc.data, nc.data.plus, along = 3)
          }
          
          
         
          env.info.yr.sub$lon.ind <- sapply(env.info.yr.sub$Lon, function(x) which.min(abs(as.numeric(dimnames(nc.data)[[1]]) - as.numeric(x))))
          env.info.yr.sub$lat.ind <- sapply(env.info.yr.sub$Lat, function(x) which.min(abs(as.numeric(dimnames(nc.data)[[2]]) - as.numeric(x))))
          env.weather <- t(sapply(seq_len(nrow(env.info.yr.sub)), function(x) nc.data[env.info.yr.sub$lon.ind[x], env.info.yr.sub$lat.ind[x], ]))
          rownames(env.weather) <- env.info.yr.sub$Environment
          env.weather <- env.weather[, 1:ndays]
          all.yrs.weather[rownames(env.weather), ] <- as.matrix(env.weather)
  
          if (verbose & sum(is.na(env.weather)) > 0) {
            NAenvs <- Envs[!stats::complete.cases(env.weather)]
            cat("\nNAs returned at ", paste(NAenvs, collapse = " "))
          }
          gc()
        }
        if (verbose) print(":)")
        suppressWarnings(file.remove(tmp.files[y]))
        return(all.yrs.weather)
    }

    if (isTRUE(ncores > 1)) { # if running in parallel
      parallel::stopCluster(cl)
      doParallel::stopImplicitCluster()
      closeAllConnections()
      if (verbose) cat("\nFinished parallel :)")
      Sys.sleep(2)
      if ("SILO_download_log.txt" %in% dir()) suppressWarnings(file.remove("SILO_download_log.txt"))
      gc(full = T)
    }
  }

  names(all.vars.weather) <- vars
  if(plus.yr){ndays <- 730}else{ndays <- 365}
  DLs <- t(sapply(Lats, function(x) {
    dls<-springpheno::daylength(daystop = 365,lat = x)
    dls<-c(dls,dls)[1:ndays]
  }))
  rownames(DLs) <- Envs
  all.vars.weather$day_length <- DLs

  NAnums <- sapply(all.vars.weather, function(x) sum(is.na(x)))
  if (verbose & sum(NAnums) > 0) {
    cat("\nNAs in:\n")
    print(NAnums)
  }

  env.info <- data.frame("Environment" = Envs, "Lat" = Lats, "Lon" = Lons)
  out <- list("data" = all.vars.weather, "Env.info" = env.info)
  return(out)
}
