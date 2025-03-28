
.onLoad <- function(libname, pkgname) {
  pages<-c("SILO"="https://s3-ap-southeast-2.amazonaws.com/silo-open-data",
           "SLGA"="https://esoil.io/TERNLandscapes/Public/Pages/SLGA/GetData-COGSDataStore.html",
           "BARRA"="https://thredds.nci.org.au/thredds/catalog/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/catalog.html",
           "cMIP6"="https://data-cbr.csiro.au/thredds/catalog/catch_all/qdc-cmip6/catalog.html")
  pges.wrkng<-sapply(pages,RCurl::url.exists)
  mssg<-paste("Connection to",names(pges.wrkng)[which(!pges.wrkng)],"API failed!")
  if(sum(!pges.wrkng)>0) packageStartupMessage(mssg)
}