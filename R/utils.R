# converts btyes to B, KB, MB or TB----
bytes.cnvrt <- function(bytes) {
  out <- ifelse(bytes < 1000, paste0(round(bytes, 2), "B"),
    ifelse(bytes < 1000000, paste0(round(bytes / 1000, 2), "KB"),
      ifelse(bytes < 1000000000, paste0(round(bytes / 1000000, 2), "MB"),
        ifelse(bytes < 1000000000000, paste0(round(bytes / 1000000000, 2), "GB"),
          paste(round(bytes / 1000000000000, 2), "TB", sep = "")
        )
      )
    )
  )
  return(out)
}


# prompt for download size----
download_data <- function(ask_before_downloading = TRUE, dl.size) {
  if (ask_before_downloading) {
    if (interactive()) {
      user_input <- readline(paste("This action will download approx ", bytes.cnvrt(dl.size), ". Would you like to proceed? (y/n): "))
      if (tolower(user_input) != "y") {
        message("Download was not performed")
        rlang::interrupt()
      }
    } else {
      message("Download was not performed, because ask_before_downloading was set to TRUE, and the environment is not interactive")
      rlang::interrupt()
    }
    message("Download beginning")
  } else {
    cat("\nTotal download approx:", bytes.cnvrt(dl.size))
    cat("\nStart time:", format(Sys.time(), "%Y-%m-%d %X"))
  }
}

# process .nc file----
nc.process <- function(nc) {
  ncfile <- ncdf4::nc_open(nc)
  nc.data <- ncdf4::ncvar_get(ncfile)
  start.day <- stringr::str_sub(string = ncfile$dim$time$units, start = nchar(ncfile$dim$time$units) - 9, end = nchar(ncfile$dim$time$units))
  dimnames(nc.data) <- list(
    ncfile$dim$lon$vals,
    ncfile$dim$lat$vals,
    as.character(as.Date(start.day) + ncfile$dim$time$vals)
  )
  ncdf4::nc_close(ncfile)
  return(nc.data)
}

# VPD function----
vpdfun <- function(tmin, tmax, relh) {
  Tave <- (tmin + tmax) / 2
  es <- 0.6108 * (exp((17.27 * Tave) / (Tave + 237.3)))
  ea <- relh / 100 * es
  VPD <- es - ea
  VPD <- VPD * 10 # Convert kpa to hpa
  return(VPD)
}

# Thermal time function----
TTfun <- function(Tci, cardT) {
  if (cardT[1] <= Tci & Tci <= cardT[2]) {
    out <- Tci
  }
  if (cardT[2] <= Tci & Tci <= cardT[3]) {
    out <- (cardT[2] / 8) * (cardT[3] - Tci)
  }
  if (Tci < cardT[1] | Tci > cardT[3]) {
    out <- 0
  }
  return(out)
}


#SMI estimation function

tmins<-weather$data$min_temp[5,]
tmaxs<-weather$data$max_temp[5,]
Tave<-(tmins+tmaxs)/2
dls<-weather$data$day_length[5,]
prs<-weather$data$daily_rain[5,]

#Potential Evapotranspiration----
PETthorn<-function(Tave,dl){
  Tm<-mean(Tave)
  N <- dl[1:length(Tave)] # day lengths
  I <- ((max(0,Tm))/5)^1.514    # heat index
  a <- (6.75e-07 * I^3) - (7.71e-05 * I^2 ) + 0.49239
  pet <- 16 * (N/360) * ((10*Tave)/I)^a
  pet[Tave < 0]<-0
  return(pet)
}

#Soil moisture index function----
SMIfun<-function(prs,tmins,tmaxs,dls,AWC,cs=2.428571,cr=0.01236147,cd=0.9644389,smi.start=0.2){
  Tave<-(tmins+tmaxs)/2
  PETs<-PETthorn(Tave = Tave,dl = dls)
  smi<-c()
  smi[1]<-smi.start
  for(i in 1:(length(prs)-1)){
    smi[i+1]<-max(0,min(1, (smi[i] / ((1+(PETs[i]/(AWC * cs))) * cd) ) + (prs[i]*cr) ))
  }
  return(smi)
}






