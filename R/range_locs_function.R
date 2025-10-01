#' @title Make a set of coordinates within the Australian cropping ranges
#'
#' @description This function uses will generate a data frame of coordinates within a given boundary of lat and lon values as well
#' as within Australian states that are within the Australian cropping ranges according to data from the Dept of Climate Change, Energy, the Environment & Water
#'
#' @param min.lon Numeric. Minimum longitude value. Default value = `110`.
#' @param max.lon Numeric. Maximum longitude value. Default value = `160`.
#' @param min.lat Numeric. Minimum latitude value. Default value = `-40`.
#' @param max.lat Numeric. Maximum latitude value. Default value = `-20`.
#' @param res Numeric. Resolution, lat/lon distance between points. Default = `0.1`.
#' @param states vector of character strings. Abbreviated names of Australian states to make points for. Options include: `c("WA","SA","NSW","QLD","VIC","TAS","NT")`
#' @returns A data frame of coordinates with Loc, Lon, and Lat column names.
#
#' @author Nick Fradgley
#'
#' @export

range.locs <- function(min.lon = 110,
                       max.lon = 160,
                       min.lat = -40,
                       max.lat = -20,
                       res = 0.1,
                       states= c("WA","SA","NSW","QLD","VIC"))
  {
  if(max.lon < min.lon) stop(crayon::red("Max lon greater than min lon"))
  if(max.lat < min.lat) stop(crayon::red("Max lat greater than min lat"))
  if(max.lon < 110) stop(crayon::red("Max lon is west of Australia"))
  if(min.lon > 160) stop(crayon::red("Min lon is east of Australia"))
  if(max.lat > -20) stop(crayon::red("Max lat is north of Australia"))
  if(min.lat < -40) stop(crayon::red("Min lat is south of Australia"))
  if(res < 0.01) warning(crayon::red("High resolution!"))
  if(res > 1) warning(crayon::red("Low resolution!"))
  if(sum(!states %in% c("WA","SA","NSW","QLD","VIC","TAS","NT"))>0) warning(crayon::red("State names are wrong!"))
  
  addrs<-"https://hub.arcgis.com/api/v3/datasets/bf09a05d02854dcd98caca1cc17b98f6_0/downloads/data?format=shp&spatialRefId=3857&where=1%3D1"
  tmp.f <- tempfile()
  tmp.f <- gsub("\\", "/", tmp.f, fixed = T)
  options(timeout = max(50000, getOption("timeout")))
  try(utils::download.file(url = addrs, destfile = tmp.f, method = "libcurl", quiet = T, mode = "wb"))
  utils::unzip(zipfile = tmp.f, exdir = gsub("\\", "/",   tempdir(), fixed = T))
  suppressWarnings(file.remove(tmp.f))
  polygons_sf <- sf::st_read(paste(gsub("\\", "/",   tempdir(), fixed = T),"/Australian_Rangeland_Boundaries.shp",sep=""),quiet=T)
  polygons_sf <- sf::st_transform(polygons_sf, crs = 4326)
  which.ranges <- which(polygons_sf$RANGELANDS==0 & polygons_sf$STATE%in%states) #Subset by states
  polygons_sf <- polygons_sf[which.ranges,] 
  
  polygons_st <- suppressWarnings(sf::st_cast(polygons_sf, "POLYGON"))
  polygon_areas <- sf::st_area(polygons_st)
  polygons_st <- polygons_st[as.numeric(polygon_areas)>1e+10,]
  
  loc.grid<-expand.grid(longitude=seq(min.lon,max.lon,res), latitude=seq(min.lat,max.lat,res))
  coordinates_df <- cbind(ID=paste(loc.grid$longitude, loc.grid$latitude, sep="_"), loc.grid) 
  coordinates_sf <- sf::st_as_sf(coordinates_df, coords = c("longitude", "latitude"), crs = 4326)
  points_in_polygons <- sf::st_join(coordinates_sf, polygons_st, join = st_intersects)
  coordinates_df <- coordinates_df[!is.na(points_in_polygons$OBJECTID),]
  colnames(coordinates_df)<-c("Loc","Lon","Lat")
  return(coordinates_df)
  }