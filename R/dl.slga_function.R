#' @title Download soil data from SLGA to a directory
#'
#' @description A function to bulk download all required TIFF files for multiple soil attributes at multiple depths
#' from the [SLGA](https://esoil.io/TERNLandscapes/Public/Pages/SLGA/index.html) data resource. TIFF files are saved in a specified directory
#' which can be used to point the [get.S.ECs()] to if used repeatedly. This will require a fast internet connection.
#'
#' @param dir Directory path that soil TIFF files will be saved in.
#' @param ncores Number (integer) of cores to use for parallel processing. Use `1` to run sequentially in series. The default (`NULL`) will
#' use the maximum available cores. If running in parallel, an output log text file will be created in the working directory.
#' @param verbose Logical. Should progress be printed? Default = TRUE.
#' @param dlprompt Logical. Should the user be prompted approve the total download size? Default = FALSE.
#'
#' @returns No value is returned but files should be saved in the `dir` directory path.
#'
#' @seealso [get.S.ECs()]
#'
#' @references
#' Grundy, M. J., Rossel, R. A. V., Searle, R. D., Wilson, P. L., Chen, C., Gregory, L. J., Grundy, M. J., Rossel,
#'      R. A. V., Searle, R. D., Wilson, P. L., Chen, C., & Gregory, L. J. (2015).
#'      [Soil and Landscape Grid of Australia. Soil Research](https://doi.org/10.1071/SR15191), 53(8), 835–844.
#'
#' @author Nick Fradgley
#'
#' @export

dl.slga<-function(dir,
                  ncores = NULL,
                  verbose = TRUE,
                  dlprompt = FALSE){
  
  is.dir<-file.info(dir)$isdir
  if(!is.dir){ stop(dir," is not a directory!") }
  
  
  atts <- unlist(c(SLGACloud::getParameterValues(Parameter = "Attribute")))[1:20]
  atts <- atts[!atts %in% c("Depth of Regolith", "Soil Organic Carbon (1\" resolution) ", "Effective Cation Exchange Capacity")]
  
  dl.size <- 700000000 * length(atts) * 6
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
    on.exit(suppressWarnings(file.remove("SLGA_soil_download_log.txt")))
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
  
  nms<-rasters$Name
  addrs<-rasters$StagingPath
  out.files<-paste(dir,"/",nms,".tif",sep="")
  
  options(timeout = max(50000, getOption("timeout")))
  
  sapply(seq_along(rasters$StagingPath),function(x){
    if (verbose == TRUE) { cat("|") }
    r <- try(terra::rast(paste("/vsicurl/",rasters$StagingPath[x], sep = "")))
    terra::writeRaster(x = r,filename = out.files[x],progress=0,overwrite=TRUE)
  })

  dl.ok<-sapply(nms,function(x) sum(as.numeric(grepl(pattern = x,x = dir(dir))))==1)
  if(sum(!dl.ok)>0){ warning("Failed to download:\n",addrs[!dl.ok])  }
  
  return(dl.ok)
  }
  if (isTRUE(ncores > 1)) { # Run in parallel
    parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
    if (verbose) {
      cat("\nFinished parallel :)\n")
    }
    Sys.sleep(2)
    closeAllConnections()
    suppressWarnings(file.remove("SLGA_soil_download_log.txt"))
  }
}
