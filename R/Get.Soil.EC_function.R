#' @title Get environmental covariates from soils data
#'
#' @description A function to derive soil Environmental Covariates (ECs) for multiple soil attributes at multiple depths
#' from the [SLGA](https://esoil.io/TERNLandscapes/Public/Pages/SLGA/index.html) data resource.
#'
#' @param Envs Vector of environment names character strings.
#' @param Lats Vector of latitude numeric values for each environment.
#' @param Lons Vector of longitude numeric values for each environment.
#' @param ncores Number (integer) of cores to use for parallel processing. Use `1` to run sequentially in series. The default (`NULL`) will
#' use the maximum available cores. If running in parallel, an output log text file will be created in the working directory.
#' @param verbose Logical. Should progress be printed?
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
                      ncores = NULL,
                      verbose = TRUE,
                      dlprompt = FALSE) {
  atts <- unlist(c(SLGACloud::getParameterValues(Parameter = "Attribute")))
  atts <- atts[1:20]
  atts <- atts[!atts %in% c("Depth of Regolith", "Soil Organic Carbon (1\" resolution) ", "Effective Cation Exchange Capacity")]
  lonlats.full <- data.frame("Loc" = paste(Lons, Lats, sep = "_"), "longitude" = Lons, "latitude" = Lats)
  lonlats.sub <- lonlats.full[!duplicated(lonlats.full$Loc), ]

  
  dl.size <- 1304 * length(atts) * 6
  download_data(dlprompt, dl.size)
  
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
      cat("\nStarting", atts[a], "\n")
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
    for (d in 1:length(depths)) {
      if (verbose == TRUE) {
        cat("\n", atts[a], "at", depths[d], "\n")
      }
      r <- NULL
      try(r <- terra::rast(paste("/vsicurl/", rasters$StagingPath[d], sep = "")))
      if (!is.null(r)) {
        vals <- unlist(terra::extract(r, as.matrix(lonlats.sub[, c("longitude", "latitude")])))
        trys <- 10
        while (sum(is.na(vals)) > 0 & trys > 0) {
          trys <- trys - 1
          na.ind <- which(is.na(vals))
          vals[is.na(vals)] <- sapply(na.ind, function(x) {
            mean(na.omit(sapply(list(
              c(0.01, 0.01),
              c(-0.01, 0.01),
              c(0.01, -0.01),
              c(-0.01, 0.01)
            ), function(c) {
              unlist(terra::extract(r, as.matrix(lonlats.sub[x, c("longitude", "latitude")]) + c))
            })))
          })
        }
        names(vals) <- lonlats.sub$Loc
        all.depth.soil[, d] <- vals[lonlats.full$Loc]
      }
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
  return(all.env.soil)
}
