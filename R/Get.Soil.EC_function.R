#' @title Get soil data from SLGA
#'
#' @description A function to derive soil Environmental Covariates (ECs) for multiple soil attributes at multiple depths
#' from the [SLGA](https://esoil.io/TERNLandscapes/Public/Pages/SLGA/index.html) data resource.
#'
#' @param Envs Vector of environment names character strings.
#' @param Lats Vector of latitude numeric values for each environment.
#' @param Lons Vector of longitude numeric values for each environment.
#' @param ncores Number (integer) of cores to use for parallel processing. Use `1` to run sequentially in series. The default (`NULL`) will
#' use the maximum available cores. If running in parallel, an output log text file will be created in the working directory.
#' @param verbose Logical. Should progress be printed? Default if TRUE.
#' @param dlprompt Logical. Should the user be prompted approve the total download size? Default it FALSE.
#'
#' @returns A data frame of soil EC values with environment names as rows and covariates as columns.
#'
#' @seealso [get.W.ECs()]
#'
#' @references
#' Grundy, M. J., Rossel, R. A. V., Searle, R. D., Wilson, P. L., Chen, C., Gregory, L. J., Grundy, M. J., Rossel,
#'      R. A. V., Searle, R. D., Wilson, P. L., Chen, C., & Gregory, L. J. (2015).
#'      [Soil and Landscape Grid of Australia. Soil Research](https://doi.org/10.1071/SR15191), 53(8), 835–844.
#'
#' @author Nick Fradgley
#'
#' @export

get.S.ECs <- function(Envs,
                      Lats,
                      Lons,
                      API.key = NULL,
                      ncores = NULL,
                      verbose = TRUE,
                      dlprompt = FALSE) {
  atts <- unlist(c(SLGACloud::getParameterValues(Parameter = "Attribute")))
  atts <- atts[1:20]
  atts <- atts[!atts %in% c("Depth of Regolith", "Soil Organic Carbon (1\" resolution) ", "Effective Cation Exchange Capacity")]
  lonlats.full <- data.frame("Loc" = paste(Lons, Lats, sep = "_"), "longitude" = Lons, "latitude" = Lats)
  lonlats.sub <- lonlats.full[!duplicated(lonlats.full$Loc), ]


  dl.size <- 300000000 * length(atts) * 6
  if (verbose) download_data(dlprompt, dl.size)

  if (is.null(ncores)) {
    ncores <- min(parallel::detectCores(), length(atts))
  }

  if (isTRUE(ncores == 1)) { # Run in series
    if (verbose) {
      cat("\nRunning in series")
    }
    `%dopar%` <- foreach::`%do%`
  }

  if (isTRUE(ncores > 1)) { # Run in parallel
    if (verbose) {
      cat("\nRunning in parallel...")
    }

    cl <- parallel::makeCluster(ncores, outfile = "SLGA_soil_download_log.txt")
    doParallel::registerDoParallel(cl)
    if (verbose) {
      cat(paste("\nProgress log output to:\n", getwd(), "/SLGA_soil_download_log.txt", sep = ""))
    }
    on.exit(closeAllConnections())
    on.exit(file.remove("SLGA_soil_download_log.txt"))
    `%dopar%` <- foreach::`%dopar%`
  }

  all.env.soil <- foreach::foreach(a = seq_along(atts), .combine = cbind, .multicombine = T) %dopar% {
    if (verbose == TRUE) {
      cat("\nStarting", atts[a])
    }
    rasters <- SLGACloud::getProductMetaData(
      Detail = "High", Attribute = atts[a],
      Component = "Modelled-Value",
      isCurrentVersion = 1
    )
    depths <- paste(rasters$UpperDepth_m, "-", rasters$LowerDepth_m, "m", sep = "")
    all.depth.soil <- matrix(NA,
      nrow = length(Envs),
      ncol = length(depths),
      dimnames = list(Envs, paste(atts[a], depths, sep = "_"))
    )

    for (d in 1:ncol(all.depth.soil)) {
      if(!is.null(API.key)){
      key <- paste0('apikey:', API.key)
      url<-gsub("https://","@",rasters$COGsPath[d])
      r <- try(terra::rast(paste("/vsicurl/https://",key,url, sep = "")),silent = T)
      if( class(r)=="try-error"){
        r <- terra::rast(paste("/vsicurl/",rasters$StagingPath[d], sep = ""))
        }
      }
      if(is.null(API.key)){
        r <- terra::rast(paste("/vsicurl/",rasters$StagingPath[d], sep = ""))
      }
      vals <- terra::extract(x = r, y = lonlats.sub[, c("longitude", "latitude")], search_radius = 100)[, 2]
      trys <- 10
      rad <- 500
      while (trys > 0 & sum(is.na(vals)) > 0) { # try filling in NAS with increasing search radius
        trys <- trys - 1
        which.nas <- which(is.na(vals))
        if (length(which.nas) < 2) {
          vals[which.nas] <- unique(terra::extract(x = r, y = lonlats.sub[which.nas, c("longitude", "latitude")], search_radius = rad)[, 2])
        }
        if (length(which.nas) > 1) {
          vals[which.nas] <- terra::extract(x = r, y = lonlats.sub[which.nas, c("longitude", "latitude")], search_radius = rad)[, 2]
        }
        rad <- rad + 500
      }
      names(vals) <- lonlats.sub$Loc
      all.depth.soil[, d] <- vals[lonlats.full$Loc]
    }
    gc()
    if (verbose & sum(!complete.cases(all.depth.soil)) > 0) {
      cat("\n NAs returned for:\n", rownames(all.depth.soil)[!complete.cases(all.depth.soil)], "\n")
    }
    return(all.depth.soil)
  }
  if (isTRUE(ncores > 1)) { # Run in parallel
    parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
    if (verbose) {
      cat("\nFinished parallel :)\n")
    }
    Sys.sleep(2)
    file.remove("SLGA_soil_download_log.txt")
    closeAllConnections()
  }

  gc(full = T)
  all.env.soil <- all.env.soil[, !colnames(all.env.soil) == "Env"]
  all.env.soil <- apply(all.env.soil, 2, function(x) {
    med <- stats::median(stats::na.omit(x))
    newx <- x
    newx[is.na(newx)] <- med
    return(newx)
  })
  
  all.env.soil <- all.env.soil[, apply(all.env.soil, 2, stats::var) > 0]
  
  isnas <- sum(is.nan(unlist(all.env.soil)) | is.na(unlist(all.env.soil)))
  if (verbose) cat(paste(isnas, "NAs returned"))
  
  if (verbose & isnas > 0) {
    cat(paste("\n NAs at:\n", paste(all.env.soil[!complete.cases(all.env.soil)], collapse = " ")))
    cat(paste("\n For:\n", paste(colnames(all.env.soil)[!complete.cases(t(all.env.soil))], collapse = " ")))
  }
  
  return(all.env.soil)
}


