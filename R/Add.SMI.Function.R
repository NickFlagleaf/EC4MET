#' @title Add estimated Soil Moisture Index (SMI) to weather data
#'
#' @description This function uses daily rainfall and temperature to estimate a SMI from a weather data object that has been output from 
#' get weather functions such as [get.SILO.weather()].
#'
#' @param weather A list of length 2. Weather data as outputted from the [get.SILO.weather()].
#'
#' `weather$data` is a list of data matrices for each covariate with rows as environments and days of the year as columns. Weather covariate names for list items should include:
#' * `daily_rain`
#' * `max_temp`
#'
#' `weather$Env.info` is a data frame of info for each environment but is not required for this function.

#' @details
#' A simplified SMI is calculated on a daily basis using a rainfall and daily maximum temperature data from the preceding 100 days. The sum of the preceding days rainfall is
#' weighted by the cumulative frying effect of the preceding days max temperatures. Values for constants have been optimised to minimise error of 
#' estimated SMI tested against more complex observed data from [TERN](https://www.tern.org.au/news-smips-soil-moisture/).
#' 
#' @returns The same `weather` object as inputed with the addition of a a matrix of daily SMI values in `$data$smi`.
#
#' @seealso [get.BARRA.weather()], [get.SILO.weather()], [get.CMIP6.weather()]
#'
#' @author Nick Fradgley
#'
#' @export

add.SMI<-function(weather){
  #Get all AWC values
  rs<-lapply(c("/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/AWC/AWC_000_005_EV_N_P_AU_TRN_N_20210614.tif",
               "/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/AWC/AWC_005_015_EV_N_P_AU_TRN_N_20210614.tif",
               "/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/AWC/AWC_015_030_EV_N_P_AU_TRN_N_20210614.tif",
               "/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/AWC/AWC_030_060_EV_N_P_AU_TRN_N_20210614.tif",
               "/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/AWC/AWC_060_100_EV_N_P_AU_TRN_N_20210614.tif",
               "/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/AWC/AWC_100_200_EV_N_P_AU_TRN_N_20210614.tif"),
             function(x) terra::rast(x) )
  AWCdata<-matrix(NA,nrow=nrow(weather$Env.info),ncol = 6,dimnames=list(weather$Env.info$Environment,c("0-0.05m","0.05-0.15m","0.15-0.3m","0.3-0.6m",
                                                                                                       "0.6-1m","1-2m")))
  for(d in 1:length(rs)){
    cat("\rGetting available water data: ",round(d/length(rs)*100,2),"%",sep="")
    vals <-  terra::extract(x = rs[[d]],y = weather$Env.info[, c("Lon", "Lat")],search_radius=5000)
    AWCdata[,d] <- vals[,2]
  }
  
  Envs<-weather$Env.info$Environment
  ndays<-100
  dayrange<-1:365
  smips<-t(sapply(Envs,function(e){
    cat("\r",rep("",50))
    cat("\rProgress: ",round(which(Envs==e)/length(Envs)*100,2),"%",sep="")
            pr<-weather$data$daily_rain[e,]
            tmax<-weather$data$max_temp[e,]
            sapply(dayrange,function(day){
              days<-(day-(ndays-1)):day
              days<-days[days>0]
              pr100<-pr[days] 
              tmax100<-tmax[days] 
              AWC<-AWCdata[e,]
              SMIfun(pr100 = pr100,tmax100 = tmax100,AWC = AWC,ndays = ndays)
              })
             }))
    colnames(smips)<-dayrange
    weather$data$smi<-smips
    return(weather)
    }



plot(weather$data$smi[10,],type="l")
barplot(weather$data$daily_rain[10,]/50,type="l",add=T,col=4,border = 4,space = 0)




test<-add.SMI(weather = wheat.area.wthr)


