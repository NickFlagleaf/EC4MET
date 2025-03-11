#' @title Get environmental covariates from soils data
#'
#' @description A function to derive soil Environmental Covariates (ECs) for multiple soil attributes at multiple depths 
#' from the [SLGA](https://esoil.io/TERNLandscapes/Public/Pages/SLGA/index.html) data resource.
#' 
#' @param Envs Vector of environment names character strings.
#' @param Lats Vector of latitude numeric values for each environment.
#' @param Lons Vector of longitude numeric values for each environment.
#' @param verbose Logical. Should progress be printed?
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
#'
#' @export

get.S.ECs <- function(Envs,
                      Lats,
                      Lons,
                      verbose = TRUE) {
  atts <- unlist(c(SLGACloud::getParameterValues(Parameter = "Attribute")))
  atts <- atts[1:20]
  atts <- atts[!atts %in% c("Depth of Regolith", "Soil Organic Carbon (1\" resolution) ", "Effective Cation Exchange Capacity")]
  lonlats <- cbind(Lons, Lats)
  colnames(lonlats) <- c("longitude", "latitude")
  all.env.soil <- data.frame("Env" = Envs)

  for (a in 1:length(atts)) {
    if (verbose == TRUE) {
      print(paste("Starting", atts[a]))
    }
    rasters <- SLGACloud::getProductMetaData(
      Detail = "High", Attribute = atts[a],
      Component = "Modelled-Value",
      isCurrentVersion = 1
    )
    depths <- paste(rasters$UpperDepth_m, "-", rasters$LowerDepth_m, "m", sep = "")
    for (d in 1:length(depths)) {
      if (verbose == TRUE) {
        print(paste("At", depths[d]))
      }
      r <- NULL
      try(r <- terra::rast(paste("/vsicurl/", rasters$StagingPath[d], sep = "")))
      if (!is.null(r)) {
        vals <- terra::extract(r, lonlats)
        colnames(vals)[1] <- paste(atts[a], depths[d], sep = "_")
        all.env.soil <- cbind.data.frame(all.env.soil, vals)
      }
      gc()
    }
  }

  rownames(all.env.soil) <- all.env.soil$Env
  all.env.soil <- all.env.soil[, !colnames(all.env.soil) == "Env"]
  all.env.soil <- apply(all.env.soil, 2, function(x) {
    med <- stats::median(stats::na.omit(x))
    newx <- x
    newx[is.na(newx)] <- med
    return(newx)
  })
  all.env.soil <- all.env.soil[, apply(all.env.soil, 2, stats::var) > 0]
  colnames(all.env.soil)[colnames(all.env.soil) %in% c(
    "Soil Organic Carbon Fractions_0-0.05m",
    "Soil Organic Carbon Fractions_0.05-0.15m",
    "Soil Organic Carbon Fractions_0.15-0.3m"
  )] <- c(
    "Soil Organic Carbon Fractions_MAOC_0-0.05m",
    "Soil Organic Carbon Fractions_MAOC0.05-0.15m",
    "Soil Organic Carbon Fractions_MAOC0.15-0.3m"
  )
  colnames(all.env.soil)[colnames(all.env.soil) %in% c(
    "Soil Organic Carbon Fractions_0-0.05m.1",
    "Soil Organic Carbon Fractions_0.05-0.15m.1"
  )] <- c(
    "Soil Organic Carbon Fractions_POC_0-0.05m",
    "Soil Organic Carbon Fractions_POC0.05-0.15m"
  )
  colnames(all.env.soil)[colnames(all.env.soil) %in% c(
    "Soil Organic Carbon Fractions_0-0.05m.2",
    "Soil Organic Carbon Fractions_0.05-0.15m.2",
    "Soil Organic Carbon Fractions_0.15-0.3m.2"
  )] <- c(
    "Soil Organic Carbon Fractions_PYOC_0-0.05m",
    "Soil Organic Carbon Fractions_PYOC_0.05-0.15m",
    "Soil Organic Carbon Fractions_PYOC_0.15-0.3m"
  )
  Smat <- all.env.soil
  return(Smat)
}
