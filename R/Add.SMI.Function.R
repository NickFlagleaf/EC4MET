#' @title Add estimated Soil Moisture Index (SMI) to weather data
#'
#' @description This function uses daily rainfall and temperature to estimate a SMI from a weather data object that has been output from
#' get weather functions such as [get.SILO.weather()].
#'
#' @param weather A list of length 2. Weather data as outputted from the [get.SILO.weather()].
#' @param verbose Logical. Should progress be printed? Default if TRUE.
#' @param API.key Optional API key for  TERN API. You can register for a key on the TERN [webpages](https://account.tern.org.au/). Default NULL
#' will download tiffs from the SLGA staging paths and may take longer to process large datasets. For large datasets with more than 100 unique locations,
#' SLGA tiff files will be temporarily downloaded rather than read directly from the API.  
#'
#' `weather$data` is a list of data matrices for each covariate with rows as environments and days of the year as columns. Weather covariate names for list items should include:
#' * `daily_rain`
#' * `max_temp`
#'
#' `weather$Env.info` is a data frame of info for each environment but is not required for this function.
#'
#' @details
#' A simplified SMI is calculated based on a dynamic water balance model where soil moisture is proportionally dried by daily evapotranspiration 
#' and increased by daily precipitation. The equation for SMI for a following day(\eqn{SMI_{i+1}}) is:
#' \deqn{SMI_{i+1} = \min{(1,\frac{SMI_i}{\frac{PET_i}{AWC \times c_s}c_d + 1} + (P_i \times c_p))}}
#' where \eqn{SMI_i} is the SMI for day \eqn{i}, \eqn{PET_i} is the potential evapotranspiration using the Thornthwaite equation for day \eqn{i}, and
#' \eqn{P_i} is the precipitation (mm) for day \eqn{i}. The constants \eqn{c_s}, \eqn{c_d}, and \eqn{c_p} were optimised based on ground truthed data from
#' [TERN](https://www.tern.org.au/news-smips-soil-moisture/)
#'
#' @returns The same `weather` list object as inputted with the addition of a matrix of daily SMI values in `...$data$smi`.
#
#' @seealso [get.BARRA.weather()], [get.SILO.weather()], [get.CMIP6.weather()]
#'
#' @author Nick Fradgley
#'
#' @export

add.SMI <- function(weather,verbose=TRUE) {
  
  rasters <- SLGACloud::getProductMetaData(
    Detail = "High", Attribute = "Available Water Capacity",
    Component = "Modelled-Value",
    isCurrentVersion = 1
  )

  Lats<-weather$Env.info$Lat
  Lons<-weather$Env.info$Lon
  lonlats.full <- data.frame("Loc" = paste(Lons, Lats, sep = "_"), "longitude" = Lons, "latitude" = Lats)
  lonlats.sub <- lonlats.full[!duplicated(lonlats.full$Loc), ]
  
  dl.n.limit<-100
  if(nrow(lonlats.sub) < dl.n.limit){
    AWCdata<-api.extrct(rasters = rasters,crds = lonlats.sub)
  }
  
  if(nrow(lonlats.sub) > dl.n.limit){
    AWCdata<-dl.extrct.tifs(addrs = rasters$StagingPath,crds = lonlats.sub)
  }
  
  colnames(AWCdata)<-paste(paste(rasters$LowerDepth_m,rasters$UpperDepth_m,sep = "-"),"m",sep="")
  AWCdata<-AWCdata[lonlats.full$Loc,]
  Envs <- weather$Env.info$Environment
  rownames(AWCdata)<-Envs
  ndays <- 100
  dayrange <- 1:365
  stps <- round(seq(1, length(Envs), length.out = 100))
  smips <- t(sapply(Envs, function(e) {
    if(verbose & e %in% Envs[stps]){
    cat("\r", rep("", 50))
    cat("\rProgress: ", round(which(Envs == e) / length(Envs) * 100), "%", sep = "")
    }
    pr <- weather$data$daily_rain[e, ]
    tmax <- weather$data$max_temp[e, ]
    tmin <- weather$data$min_temp[e, ]
    dl <- weather$data$day_length[e, ]
    awc <- mean(AWCdata[e, ])
    smi <- SMIfun(prs = pr, tmins = tmin, tmaxs = tmax, dls = dl, AWC = awc)
    return(smi)
  }))
  colnames(smips) <- dayrange
  weather$data$smi <- smips
  return(weather)
}
