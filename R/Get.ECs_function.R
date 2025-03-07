#' @title Get environmental covariates from weather data
#'
#' @description A function to calculate Environmental Covariates (ECs) from daily weather data such as derived from the [get.SILO.weather()] function.
#'
#' @param weather A two level list of Weather data as outputted from the [get.SILO.weather()]. `Weather$data` is a list of data matrices for each covariate
#' with rows as environments and days of the year as columns. Weather covariate names should be:
#' 
#' *"daily_rain"
#' 
#' *"max_temp"
#' 
#' *"min_temp"
#' 
#' *"vp_deficit"
#' 
#' *"radiation"
#'
#' `obs.wthr$Env.info` is a data frame of info for each environment that includes a `Lat` columnof latitude values for which day lengths are calculated.
#' @param  sow.dates Vector of character strings of dates of sowing for each trail environment in dd/mm/yyy format. Must be in the same order as the
#' rownames of weather data in `weather$data` matrices.
#' @param cardT Optional. Minimum, optimal and maximum cardinal temperatures to calculate thermal time. Default values are min = 0, opt = 26, and max = 34.
#' Custom values can be used to define other crop phenologies and growth rates
#' @param stg.TT Optional. Estimated thermal time parameters between wheat crop growth stages. Default values:
#'
#'    * Emergence to End of juvenile growth stages = 500
#'
#'    * Heading to Flowering = 250
#'
#'    * Flowering to Start of grain fill = 250
#'
#'    * Start of Grain fill to End of grain fill = 250
#'
#'    * End of grain fill to Maturity = 400
#'
#'    Custom values can be used to define other crop phenologies and growth rates
#'
#' @param DTH.TT Optional. Estimated thermal time from sowing that heading growth stage occurs. Default value is 1285.
#' @param verbose Logical. Should progress be printed?
#'
#' @returns A data frame of weather EC values with environment names as rows and covariates as columns.
#' @export

get.W.ECs <- function(weather,
                      sow.dates,
                      cardT = c(0, 26, 34),
                      stg.TT = c(500, 250, 250, 250, 400),
                      DTH.TT = 1285,
                      verbose = TRUE) {
  TTfun <- function(Tci) {
    if (0 <= Tci & Tci <= 26) {
      out <- Tci
    }
    if (26 <= Tci & Tci <= 34) {
      out <- (26 / 8) * (34 - Tci)
    }
    if (Tci < 0 | Tci > 34) {
      out <- 0
    }
    return(out)
  }

  Tc <- (weather$data$max_temp + weather$data$min_temp) / 2
  all.envDailyTT <- t(apply(Tc, 1, function(x) sapply(x, TTfun)))

  stage.names <- c("Sowing", "Emergence", "End of juvenile", "Heading", "Flowering", "Start of grain filling", "End of grain filling", "Maturity")
  Envs <- rownames(weather$data$daily_rain)
  all.env.stages <- matrix(NA,
    nrow = length(Envs), ncol = length(stage.names),
    dimnames = list(Envs, stage.names)
  )
  sow.dates <- as.Date(sow.dates, format = "%d/%m/%y")
  yrday1 <- as.Date(paste(stringr::str_sub(string = sow.dates, 1, 4), "-01-01", sep = ""))
  sowdays <- sow.dates - yrday1
  cat("\nStarting growth stage estimates\n")
  for (i in 1:length(Envs)) {
    DailyTT <- all.envDailyTT[i, ]
    Dailypr <- weather$data$daily_rain[i, ]
    sow.day <- sowdays[i]
    stages <- rep(NA, 8)
    names(stages) <- stage.names
    stages[1] <- 1
    wetsow <- sum(Dailypr[c(sow.day - 7):sow.day]) > 3
    if (wetsow) {
      stages[2] <- 14
    } else {
      stages[2] <- min(min(which(Dailypr[sow.day:364] > 0)) + 14, 40)
    }
    stages[3] <- which.min(abs(cumsum(DailyTT[sow.day:364]) - (cumsum(DailyTT[sow.day:364])[stages[2]] + stg.TT[1])))
    stages[4] <- which.min(abs(cumsum(DailyTT[sow.day:364]) - DTH.TT))
    stages[5] <- which.min(abs(cumsum(DailyTT[sow.day:364]) - (cumsum(DailyTT[sow.day:364])[stages[4]] + stg.TT[2])))
    stages[6] <- which.min(abs(cumsum(DailyTT[sow.day:364]) - (cumsum(DailyTT[sow.day:364])[stages[5]] + stg.TT[3])))
    stages[7] <- which.min(abs(cumsum(DailyTT[sow.day:364]) - (cumsum(DailyTT[sow.day:364])[stages[6]] + stg.TT[4])))
    stages[8] <- which.min(abs(cumsum(DailyTT[sow.day:364]) - (cumsum(DailyTT[sow.day:364])[stages[7]] + stg.TT[5])))
    if (verbose == TRUE & i %in% round(seq(1, length(Envs), length.out = 100))) {
      cat("|", round(i / length(Envs) * 100), "%", sep = "")
    }
    all.env.stages[i, ] <- stages
  }
  all.env.stages <- as.data.frame(all.env.stages)


  # Define stress covariates----
  interval.names <- c("Sow2Emer", "Emer2Juv", "Juv2He", "He2Flw", "Flw2Sgf", "Sgf2Egf", "Egf2mat")

  # N days per stage----
  if (verbose == TRUE) {
    cat("\nStarting N days per stage")
  }
  Ndays.per.stage <- t(sapply(1:length(Envs), function(e) sapply(2:ncol(all.env.stages), function(s) all.env.stages[e, s] - all.env.stages[e, s - 1])))
  colnames(Ndays.per.stage) <- paste("Ndays_", interval.names, sep = "")

  # N days from sowing to flowering----
  if (verbose == TRUE) {
    cat("\nStarting N days from sowing to flowering")
  }
  Ndays_Sow2Flw <- all.env.stages[, "Flowering"]
  names(Ndays_Sow2Flw) <- Envs

  # N days from flowering to end of grain fill----
  if (verbose == TRUE) {
    cat("\nStarting N days from flowering to end of grain fill")
  }
  Ndays_Flw2Egf <- all.env.stages[, "End of grain filling"] - all.env.stages[, "Flowering"]
  names(Ndays_Flw2Egf) <- Envs

  # total rain per stage-----
  if (verbose == TRUE) {
    cat("\nStarting total rain per stage\n")
  }
  {
    Sum.rain.per.stage <- matrix(NA,
      nrow = length(Envs), ncol = length(interval.names),
      dimnames = list(Envs, paste("TotRain_", interval.names, sep = ""))
    )
    e <- 1
    stps <- round(seq(1, length(Envs), length.out = 100))
    for (e in e:length(Envs)) {
      yr.pr <- unlist(weather$data$daily_rain[e, ])
      stgs <- unlist(all.env.stages[e, ])
      Sum.rain.per.stage[e, ] <- sapply(2:length(stgs), function(s) {
        sum(yr.pr[(sow.day + stgs[s - 1]):(sow.day + stgs[s])])
      })
      if (verbose == TRUE & e %in% stps) {
        cat("|", round(e / length(Envs) * 100), "%", sep = "")
      }
    }
  }


  # total rain from sowing to flowering---------
  if (verbose == TRUE) {
    cat("\nStarting total rain from sowing to flowering")
  }
  TotRain_Sow2Flw <- sapply(1:length(Envs), function(e) {
    yr.pr <- unlist(weather$data$daily_rain[e, ])
    sum.pr <- sum(yr.pr[sow.day:364][all.env.stages$Sowing[e]:all.env.stages$Flowering[e]])
    return(sum.pr)
  })
  names(TotRain_Sow2Flw) <- Envs

  # total rain from flowering to end of grain fill--------
  if (verbose == TRUE) {
    cat("\nStarting total rain from flowering to end of grain fill")
  }
  TotRain_Flw2Egf <- sapply(1:length(Envs), function(e) {
    yr.pr <- unlist(weather$data$daily_rain[e, ])
    sum.pr <- sum(yr.pr[sow.day:364][all.env.stages$Flowering[e]:all.env.stages$`End of grain filling`[e]])
    return(sum.pr)
  })
  names(TotRain_Flw2Egf) <- Envs

  # Total stored soil moisture before sowing---------
  if (verbose == TRUE) {
    cat("\nStarting Total stored soil moisture before sowing")
  }
  TotRain_priorSow <- sapply(1:length(Envs), function(e) {
    yr.pr <- unlist(weather$data$daily_rain[e, ])
    sum.pr <- sum(yr.pr[1:sow.day])
    return(sum.pr)
  })
  names(TotRain_priorSow) <- Envs


  # mean temp per stage------
  if (verbose == TRUE) {
    cat("\nStarting mean temp per stage\n")
  }
  mean.temp.per.stage <- matrix(NA,
    nrow = length(Envs), ncol = length(interval.names),
    dimnames = list(Envs, paste("Avtemp_", interval.names, sep = ""))
  )
  stps <- round(seq(1, length(Envs), length.out = 100))
  for (e in 1:length(Envs)) {
    yr.tmax <- unlist(weather$data$max_temp[e, ])
    yr.tmin <- unlist(weather$data$min_temp[e, ])
    stgs <- unlist(all.env.stages[e, ])
    mean.temp.per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) mean(((yr.tmax[(sow.day + stgs[s - 1]):(sow.day + stgs[s])] + yr.tmin[(sow.day + stgs[s - 1]):(sow.day + stgs[s])]) / 2)))
    if (verbose == TRUE & e %in% stps) {
      cat("|", round(e / length(Envs) * 100), "%", sep = "")
    }
  }

  # mean max temp per stage------
  if (verbose == TRUE) {
    cat("\nStarting mean max temp per stage\n")
  }
  mean.max.temp.per.stage <- matrix(NA,
    nrow = length(Envs), ncol = length(interval.names),
    dimnames = list(Envs, paste("Avmaxtemp_", interval.names, sep = ""))
  )
  stps <- round(seq(1, length(Envs), length.out = 100))
  for (e in 1:length(Envs)) {
    stgs <- unlist(all.env.stages[e, ])
    yr.tmax <- unlist(weather$data$max_temp[e, ])
    mean.max.temp.per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) mean(yr.tmax[(sow.day + stgs[s - 1]):(sow.day + stgs[s])]))
    if (verbose == TRUE & e %in% stps) {
      cat("|", round(e / length(Envs) * 100), "%", sep = "")
    }
  }

  # mean min temp per stage-----
  if (verbose == TRUE) {
    cat("\nStarting mean min temp per stage\n")
  }
  mean.min.temp.per.stage <- matrix(NA,
    nrow = length(Envs), ncol = length(interval.names),
    dimnames = list(Envs, paste("AvMintemp_", interval.names, sep = ""))
  )
  stps <- round(seq(1, length(Envs), length.out = 100))
  for (e in 1:length(Envs)) {
    stgs <- unlist(all.env.stages[e, ])
    yr.tmin <- unlist(weather$data$min_temp[e, ])
    mean.min.temp.per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) mean(yr.tmin[(sow.day + stgs[s - 1]):(sow.day + stgs[s])]))
    if (verbose == TRUE & e %in% stps) {
      cat("|", round(e / length(Envs) * 100), "%", sep = "")
    }
  }


  # Sum dry days per stage-------
  if (verbose == TRUE) {
    cat("\nStarting Sum dry days per stage\n")
  }
  Sum.drydays.per.stage <- matrix(NA,
    nrow = length(Envs), ncol = length(interval.names),
    dimnames = list(Envs, paste("Ndd_", interval.names, sep = ""))
  )
  stps <- round(seq(1, length(Envs), length.out = 100))
  for (e in 1:length(Envs)) {
    stgs <- unlist(all.env.stages[e, ])
    yr.pr <- unlist(weather$data$daily_rain[e, ])
    Sum.drydays.per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) sum(yr.pr[(sow.day + stgs[s - 1]):(sow.day + stgs[s])] < 1))
    if (verbose == TRUE & e %in% stps) {
      cat("|", round(e / length(Envs) * 100), "%", sep = "")
    }
  }

  # Sum minT < 0C per stage---------
  if (verbose == TRUE) {
    cat("\nStarting Sum minT < 0C per stage\n")
  }
  MinTbelow0per.stage <- matrix(NA,
    nrow = length(Envs), ncol = length(interval.names),
    dimnames = list(Envs, paste("Ndays<0_", interval.names, sep = ""))
  )
  stps <- round(seq(1, length(Envs), length.out = 100))
  for (e in 1:length(Envs)) {
    stgs <- unlist(all.env.stages[e, ])
    yr.tmin <- unlist(weather$data$min_temp[e, ])
    MinTbelow0per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) sum(yr.tmin[(sow.day + stgs[s - 1]):(sow.day + stgs[s])] < 0))
    if (verbose == TRUE & e %in% stps) {
      cat("|", round(e / length(Envs) * 100), "%", sep = "")
    }
  }


  # Frosts at flowering-----
  if (verbose == TRUE) {
    cat("\nStarting Frosts at flowering")
  }
  MinTbelow0FLW <- c()
  for (e in 1:length(Envs)) {
    yr.tmin <- unlist(weather$data$min_temp[e, ])
    MinTbelow0FLW[e] <- sum(c(yr.tmin[sow.day:364] < 0)[(all.env.stages$Flowering[e] - 7):(all.env.stages$Flowering[e] + 7)])
  }
  names(MinTbelow0FLW) <- Envs


  # warm days over 26c------
  if (verbose == TRUE) {
    cat("\nStarting warm days over 26c\n")
  }
  MaxToverr26per.stage <- matrix(NA,
    nrow = length(Envs), ncol = length(interval.names),
    dimnames = list(Envs, paste("Ndays>26_", interval.names, sep = ""))
  )
  stps <- round(seq(1, length(Envs), length.out = 100))
  for (e in 1:length(Envs)) {
    stgs <- unlist(all.env.stages[e, ])
    yr.tmax <- unlist(weather$data$max_temp[e, ])
    MaxToverr26per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) sum(yr.tmax[(sow.day + stgs[s - 1]):(sow.day + stgs[s])] > 26))
    if (verbose == TRUE & e %in% stps) {
      cat("|", round(e / length(Envs) * 100), "%", sep = "")
    }
  }


  # Hot days over 34 per stage-------
  if (verbose == TRUE) {
    cat("\nStarting Hot days per stage\n")
  }
  MaxToverr34per.stage <- matrix(NA,
    nrow = length(Envs), ncol = length(interval.names),
    dimnames = list(Envs, paste("Ndays>34_", interval.names, sep = ""))
  )
  stps <- round(seq(1, length(Envs), length.out = 100))
  for (e in 1:length(Envs)) {
    stgs <- unlist(all.env.stages[e, ])
    yr.tmax <- unlist(weather$data$max_temp[e, ])
    MaxToverr34per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) sum(yr.tmax[(sow.day + stgs[s - 1]):(sow.day + stgs[s])] > 34))
    if (verbose == TRUE & e %in% stps) {
      cat("|", round(e / length(Envs) * 100), "%", sep = "")
    }
  }

  # mean sunshine per stage-------
  if (verbose == TRUE) {
    cat("\nStarting SR days per stage\n")
  }
  AveSR.per.stage <- matrix(NA,
    nrow = length(Envs), ncol = length(interval.names),
    dimnames = list(Envs, paste("AveSR_", interval.names, sep = ""))
  )
  stps <- round(seq(1, length(Envs), length.out = 100))
  for (e in 1:length(Envs)) {
    stgs <- unlist(all.env.stages[e, ])
    yr.sr <- unlist(weather$data$radiation[e, ])
    AveSR.per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) sum(yr.sr[(sow.day + stgs[s - 1]):(sow.day + stgs[s])] > 34))
    if (verbose == TRUE & e %in% stps) {
      cat("|", round(e / length(Envs) * 100), "%", sep = "")
    }
  }

  # mean VPD per stage-------
  if (verbose == TRUE) {
    cat("\nStarting SR days per stage\n")
  }
  AveVPD.per.stage <- matrix(NA,
    nrow = length(Envs), ncol = length(interval.names),
    dimnames = list(Envs, paste("AveVPD_", interval.names, sep = ""))
  )
  stps <- round(seq(1, length(Envs), length.out = 100))
  for (e in 1:length(Envs)) {
    stgs <- unlist(all.env.stages[e, ])
    yr.vpd <- unlist(weather$data$vp_deficit[e, ])
    AveVPD.per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) mean(yr.vpd[(sow.day + stgs[s - 1]):(sow.day + stgs[s])]))
    if (verbose == TRUE & e %in% stps) {
      cat("|", round(e / length(Envs) * 100), "%", sep = "")
    }
  }

  # mean day lengths per stage---------
  if (verbose == TRUE) {
    cat("\nStarting mean day lengths per stage\n")
  }
  all.lats <- weather$Env.info$Lat
  MmeanDLper.stage <- matrix(NA,
    nrow = length(Envs), ncol = length(interval.names),
    dimnames = list(Envs, paste("AveDL_", interval.names, sep = ""))
  )
  stps <- round(seq(1, length(Envs), length.out = 100))
  for (e in 1:length(Envs)) {
    lat <- all.lats[e]
    yr.DLs <- chillR::daylength(latitude = lat, JDay = 1:364, notimes.as.na = FALSE)$Daylength
    stgs <- unlist(all.env.stages[e, ])
    MmeanDLper.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) mean(yr.DLs[as.numeric(sow.day + stgs[(s - 1):s])]))
    if (verbose == TRUE & e %in% stps) {
      cat("|", round(e / length(Envs) * 100), "%", sep = "")
    }
  }

  # Make weather matrix------
  Wmat <- cbind.data.frame(Ndays.per.stage,
    Ndays_Sow2Flw,
    Ndays_Flw2Egf,
    Sum.rain.per.stage,
    TotRain_Sow2Flw,
    TotRain_Flw2Egf,
    TotRain_priorSow,
    mean.temp.per.stage,
    mean.min.temp.per.stage,
    "Mintemp<0_Flw" = MinTbelow0FLW,
    mean.max.temp.per.stage,
    Sum.drydays.per.stage,
    MinTbelow0per.stage,
    MaxToverr26per.stage,
    MaxToverr34per.stage,
    AveSR.per.stage,
    AveVPD.per.stage,
    MmeanDLper.stage
  )
  rownames(Wmat) <- Envs

  if (verbose == TRUE) {
    cat(paste("\n", sum(is.nan(unlist(Wmat)) | is.na(unlist(Wmat))), "missing values"))
  }

  return(Wmat)
}
