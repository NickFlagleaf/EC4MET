#' @title Get environmental covariates from weather data
#'
#' @description A function to calculate Environmental Covariates (ECs) from daily weather data such as derived from the [get.SILO.weather()] or [get.BARRA.weather()] functions.
#'
#' @param weather A list of length 2. Weather data as outputted from the [get.SILO.weather()].
#'
#' `weather$data` is a list of data matrices for each covariate with rows as environments and days of the year as columns. Weather covariate names for list items should be:
#' * `daily_rain`
#' * `max_temp`
#' * `min_temp`
#' * `vp_deficit`
#' * `radiation`
#'
#' `weather$Env.info` is a data frame of info for each environment that includes a `Lat` column of latitude values for which day lengths are calculated.
#'
#' @param  sow.dates Vector of character strings of dates of sowing for each trail environment in dd/mm/yyy format. Must be in the same order as the
#' rownames of weather data in `weather$data` matrices.
#' @param cardT Optional. Vector of minimum, optimal and maximum cardinal temperatures to calculate thermal time. Default values are min = 0, opt = 26, and max = 34.
#' Custom values can be used to define other crop growth parameters.
#' @param stg.TT Optional. Thermal time parameters that are used to estimate the number of days between wheat crop growth stages (see details). Default values:
#' * Emergence to End of juvenile growth stages = 500
#' * Heading to Flowering = 250
#' * Flowering to Start of grain fill = 250
#' * Start of Grain fill to End of grain fill = 250
#' * End of grain fill to Maturity = 400
#'
#' Custom values can be used to define other crop phenologies and growth rates
#'
#' @param DTH.TT Optional. Estimated thermal time from sowing that heading growth stage occurs. Default value is 1285.
#' @param verbose Logical. Should progress be printed?
#'
#' @details
#' ECs are calculated for periods between crop growth stages that are estimated based on a thermal time degree days model defined by the `cardT` parameters.
#' Crop growth stages abbreviations and equivalent Zadocks scale:
#' * `Sow` - Sowing (GS0)
#' * `Emer` - Emergence (GS10)
#' * `Juv` - End of Juvenile (GS30)
#' * `He` - Heading (GS55)
#' * `Flow` - Flowering (GS65)
#' * `Sgf` - Start of grain filling (GS71)
#' * `Egf` - End of Grain filling (GS87)
#' * `Mat` - Maturity (GS92)
#'
#' Other abbreviations for ECs calculated between growth stage intervals and combined into EC names include:
#' * `Ndays` - Number of days
#' * `TotRain` - Total rainfall (mm)
#' * `Avtemp` - Average of average daily temperatures (°C)
#' * `AvMintemp` - Average of daily minimum temperatures (°C)
#' * `Avmaxtemp` -  Average of daily minimum temperatures (°C)
#' * `Ndays<0` - Number of frost days that the minimum temperature was below 0 °C
#' * `Ndays>26` - Number of warm days that max temp was over 26 °C
#' * `Ndays>34` - Number of hot days that max temp was over 34 °C
#' * `AveSR` -  Average solar radiation (MJ m<sup>-2</sup>)
#' * `AveVPD` - Average vapour pressure deficit (hPa)
#' * `AvePQ` - Average photothermal quotient (MJ m<sup>-2</sup> day<sup>-1</sup> °C<sup>-1</sup>)
#' * `AveDL` - Average day length (hr)
#'
#' For example, `TotRain_He2Flw` indicates the total rainfall between the heading and flowering growth stages.
#'
#' Other specific ECs include:
#' * `TotRain_priorSow` - The total rainfall between Jan 1<sup>st</sup> and the sowing day (mm)
#' * `Mintemp<0_Flw` - Number number of frost days within 7 days of the estimated flowering date
#'
#' For details of how ECs are calculated, see Fradgley et al. 2025.
#'
#' @returns A list of length 2:
#' * `$ECs` - A data frame of weather EC values with environment names as rows and covariates as columns.
#' * `$gs.dates` - A data frame of estimated dates in yyyy-mm-dd format of each growth stage per environment with environment names as rows and abbreviated
#' growth stage names as columns.
#'
#' @seealso [get.S.ECs()], [get.BARRA.weather()], [get.SILO.weather()]
#'
#' @references
#' * Fradgley et al. (2025) Prediction of Australian wheat genotype by environment interactions and mega-environments,
#'  Under review.
#' * Zadoks, J. C., Chang, T. T., & Konzak, C. F. (1974). [A decimal code for the growth stages of cereals](https://doi.org/10.1111/j.1365-3180.1974.tb01084.x).
#'    Weed research, 14(6), 415-421.
#'
#' @export

get.W.ECs <- function(weather,
                      sow.dates,
                      cardT = NULL,
                      stg.TT = NULL,
                      DTH.TT = NULL,
                      verbose = TRUE) {
  if (is.null(cardT)) {
    cardT <- c("min" = 0, "opt" = 26, "max" = 34)
  }

  if (is.null(DTH.TT)) {
    DTH.TT <- 1285
  }

  if (is.null(stg.TT)) {
    stg.TT <- c(
      "Sow-Juv" = 500,
      "He-Flow" = 250,
      "Flow-Sgf" = 250,
      "Sgf-Egf" = 250,
      "Egf-Mat" = 400
    )
  }

  TTfun <- function(Tci) {
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


  Tc <- (weather$data$max_temp + weather$data$min_temp) / 2
  all.envDailyTT <- t(apply(Tc, 1, function(x) sapply(x, TTfun)))

  stage.names <- c("Sow", "Emer", "Juv", "He", "Flow", "Sgf", "Egf", "Mat")
  Envs <- rownames(weather$data$daily_rain)
  all.env.stages <- matrix(NA,
    nrow = length(Envs), ncol = length(stage.names),
    dimnames = list(Envs, stage.names)
  )
  sow.dates <- as.Date(sow.dates, tryFormats = c("%d/%m/%Y", "%Y/%m/%d", "%d-%m-%Y", "%Y-%m-%d"))
  yrday1 <- as.Date(paste(stringr::str_sub(string = sow.dates, 1, 4), "-01-01", sep = ""))
  sowdays <- sow.dates - yrday1

  if (verbose & !length(sow.dates) == length(Envs)) {
    cat("\n!length(sow.dates)==length(Envs)")
  }

  if (verbose & sum(is.na(sow.dates)) > 0) {
    cat(paste("\nError for sow.dates at ", paste(Envs[is.na(sow.dates)], collapse = " "), sep = ""))
  }

  if (verbose) {
    cat("\nStarting growth stage estimates ")
  }

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

  {
    if (verbose == TRUE) {
      cat("\nStarting days per stage ")
    }
    # N days per stage----
    Ndays.per.stage <- sapply(2:ncol(all.env.stages), function(s) all.env.stages[, s] - all.env.stages[, s - 1])
    colnames(Ndays.per.stage) <- paste("Ndays_", interval.names, sep = "")

    # N days from sowing to flowering----
    Ndays_Sow2Flw <- all.env.stages[, "Flow"]

    # N days from flowering to end of grain fill----
    Ndays_Flw2Egf <- all.env.stages[, "Egf"] - all.env.stages[, "Flow"]

    ndays.ECs <- cbind(Ndays.per.stage, Ndays_Sow2Flw, Ndays_Flw2Egf)
  }

  { # total rain per stage-----
    if (verbose == TRUE) {
      cat("\nStarting total and rain and dry days per stages ")
    }
    Sum.rain.per.stage <- matrix(NA,
      nrow = length(Envs), ncol = length(interval.names) + 3,
      dimnames = list(Envs, c(
        paste("TotRain_", interval.names, sep = ""),
        "TotRain_Sow2Flw", "TotRain_Flw2Egf", "TotRain_priorSow"
      ))
    )
    Sum.drydays.per.stage <- matrix(NA,
      nrow = length(Envs), ncol = length(interval.names),
      dimnames = list(Envs, paste("Ndd_", interval.names, sep = ""))
    )
    e <- 1
    stps <- round(seq(1, length(Envs), length.out = 100))
    for (e in e:length(Envs)) {
      yr.pr <- unlist(weather$data$daily_rain[e, ])
      stgs <- unlist(all.env.stages[e, ])
      sowday <- sowdays[e]
      Sum.rain.per.stage[e, paste("TotRain_", interval.names, sep = "")] <- sapply(2:length(stgs), function(s) {
        sum(yr.pr[(sowday + stgs[s - 1]):(sowday + stgs[s])])
      })

      Sum.rain.per.stage[e, "TotRain_Sow2Flw"] <- sum(yr.pr[sowday:(sowday + all.env.stages$Flow[e])])
      Sum.rain.per.stage[e, "TotRain_Flw2Egf"] <- sum(yr.pr[(sowday + all.env.stages$Flow[e]):(sowday + all.env.stages$Egf[e])])
      Sum.rain.per.stage[e, "TotRain_priorSow"] <- sum(yr.pr[1:sowday])

      Sum.drydays.per.stage[e, ] <- sapply(2:length(stgs), function(s) sum(yr.pr[(sowday + stgs[s - 1]):(sowday + stgs[s])] < 1))

      if (verbose == TRUE & e %in% stps) {
        cat("|", round(e / length(Envs) * 100), "%", sep = "")
      }
    }
    rain.ECs <- cbind(Sum.rain.per.stage, Sum.drydays.per.stage)
  }

  { # mean temp per stage------
    if (verbose == TRUE) {
      cat("\nStarting temps per stage ")
    }
    mean.temp.per.stage <- matrix(NA,
      nrow = length(Envs), ncol = length(interval.names),
      dimnames = list(Envs, paste("Avtemp_", interval.names, sep = ""))
    )
    mean.max.temp.per.stage <- matrix(NA,
      nrow = length(Envs), ncol = length(interval.names),
      dimnames = list(Envs, paste("Avmaxtemp_", interval.names, sep = ""))
    )
    mean.min.temp.per.stage <- matrix(NA,
      nrow = length(Envs), ncol = length(interval.names),
      dimnames = list(Envs, paste("AvMintemp_", interval.names, sep = ""))
    )
    MinTbelow0per.stage <- matrix(NA,
      nrow = length(Envs), ncol = length(interval.names) + 1,
      dimnames = list(Envs, c(paste("Ndays<0_", interval.names, sep = ""), "Ndays<0_Flw"))
    )
    MaxToverr26per.stage <- matrix(NA,
      nrow = length(Envs), ncol = length(interval.names),
      dimnames = list(Envs, paste("Ndays>26_", interval.names, sep = ""))
    )
    MaxToverr34per.stage <- matrix(NA,
      nrow = length(Envs), ncol = length(interval.names),
      dimnames = list(Envs, paste("Ndays>34_", interval.names, sep = ""))
    )


    mean.temp <- (weather$data$max_temp + weather$data$min_temp) / 2
    stps <- round(seq(1, length(Envs), length.out = 100))
    for (e in 1:length(Envs)) {
      yr.tmean <- unlist(mean.temp[e, ])
      yr.tmax <- unlist(weather$data$max_temp[e, ])
      yr.tmin <- unlist(weather$data$min_temp[e, ])
      stgs <- unlist(all.env.stages[e, ])
      sowday <- sowdays[e]
      mean.temp.per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) mean(mean.temp[(sowday + stgs[s - 1]):(sowday + stgs[s])]))
      mean.max.temp.per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) mean(yr.tmax[(sowday + stgs[s - 1]):(sowday + stgs[s])]))
      mean.min.temp.per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) mean(yr.tmin[(sowday + stgs[s - 1]):(sowday + stgs[s])]))
      MinTbelow0per.stage[e, paste("Ndays<0_", interval.names, sep = "")] <- sapply(
        2:ncol(all.env.stages),
        function(s) sum(yr.tmin[(sowday + stgs[s - 1]):(sowday + stgs[s])] < 0)
      )
      MinTbelow0per.stage[e, "Ndays<0_Flw"] <- sum(c(yr.tmin[sowday:364] < 0)[(all.env.stages$Flow[e] - 7):(all.env.stages$Flow[e] + 7)])
      MaxToverr26per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) sum(yr.tmax[(sowday + stgs[s - 1]):(sowday + stgs[s])] > 26))
      MaxToverr34per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) sum(yr.tmax[(sowday + stgs[s - 1]):(sowday + stgs[s])] > 34))

      if (verbose == TRUE & e %in% stps) {
        cat("|", round(e / length(Envs) * 100), "%", sep = "")
      }
    }
    temp.ECS <- cbind(
      mean.temp.per.stage, mean.max.temp.per.stage, mean.min.temp.per.stage,
      MinTbelow0per.stage, MaxToverr26per.stage, MaxToverr34per.stage
    )
  }

  { # mean sunshine per stage-------
    if (verbose == TRUE) {
      cat("\nStarting sol rad and PQ per stage ")
    }
    AveSR.per.stage <- matrix(NA,
      nrow = length(Envs), ncol = length(interval.names),
      dimnames = list(Envs, paste("AveSR_", interval.names, sep = ""))
    )
    AvePQ.per.stage <- matrix(NA,
      nrow = length(Envs), ncol = length(interval.names),
      dimnames = list(Envs, paste("AvePQ_", interval.names, sep = ""))
    )
    stps <- round(seq(1, length(Envs), length.out = 100))
    for (e in 1:length(Envs)) {
      stgs <- unlist(all.env.stages[e, ])
      yr.sr <- unlist(weather$data$radiation[e, ])
      yr.tmean <- unlist(mean.temp[e, ])
      sowday <- sowdays[e]
      AveSR.per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) mean(yr.sr[(sowday + stgs[s - 1]):(sowday + stgs[s])]))
      AvePQ.per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) {
        mean((yr.sr[(sowdays[e] + stgs[s - 1]):(sowdays[e] + stgs[s])] * .47) / yr.tmean[(sowdays[e] + stgs[s - 1]):(sowdays[e] + stgs[s])])
      })

      if (verbose == TRUE & e %in% stps) {
        cat("|", round(e / length(Envs) * 100), "%", sep = "")
      }
    }
    sunECs <- cbind(AveSR.per.stage, AvePQ.per.stage)
  }

  { # mean VPD per stage-------
    if (verbose == TRUE) {
      cat("\nStarting VPD days per stage ")
    }
    AveVPD.per.stage <- matrix(NA,
      nrow = length(Envs), ncol = length(interval.names),
      dimnames = list(Envs, paste("AveVPD_", interval.names, sep = ""))
    )

    stps <- round(seq(1, length(Envs), length.out = 100))
    for (e in 1:length(Envs)) {
      stgs <- unlist(all.env.stages[e, ])
      yr.vpd <- unlist(weather$data$vp_deficit[e, ])
      AveVPD.per.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) mean(yr.vpd[(sowdays[e] + stgs[s - 1]):(sowdays[e] + stgs[s])]))
      if (verbose == TRUE & e %in% stps) {
        cat("|", round(e / length(Envs) * 100), "%", sep = "")
      }
    }
  }

  # mean day lengths per stage---------
  if (verbose == TRUE) {
    cat("\nStarting mean day lengths per stage ")
  }
  MeanDLper.stage <- matrix(NA,
    nrow = length(Envs), ncol = length(interval.names),
    dimnames = list(Envs, paste("AveDL_", interval.names, sep = ""))
  )
  stps <- round(seq(1, length(Envs), length.out = 100))
  for (e in 1:length(Envs)) {
    yr.DLs <- unlist(weather$data$day_length[e, ])
    stgs <- unlist(all.env.stages[e, ])
    MeanDLper.stage[e, ] <- sapply(2:ncol(all.env.stages), function(s) mean(yr.DLs[as.numeric(sowdays[e] + stgs[(s - 1):s])]))
    if (verbose == TRUE & e %in% stps) {
      cat("|", round(e / length(Envs) * 100), "%", sep = "")
    }
  }



  # Make weather matrix------
  Wmat <- cbind.data.frame(ndays.ECs, rain.ECs, temp.ECS, sunECs, AveVPD.per.stage, MeanDLper.stage)

  rownames(Wmat) <- Envs

  isnas <- sum(is.nan(unlist(Wmat)) | is.na(unlist(Wmat)))
  if (verbose) {
    cat(paste("\n", isnas, "NAs"))
  }
  if (verbose & isnas > 0) {
    cat(paste("\n NAs at:\n", paste(Envs[!complete.cases(Wmat)], collapse = " ")))
    cat(paste("\n For:\n", paste(colnames(Wmat)[!complete.cases(t(Wmat))], collapse = " ")))
  }


  GS.dates <- t(sapply(1:nrow(all.env.stages), function(x) as.character(sow.dates[x] + unlist(all.env.stages[x, ]) - 1)))
  dimnames(GS.dates) <- list(Envs, stage.names)

  out <- list(
    "gs.dates" = GS.dates,
    "ECs" = Wmat
  )

  return(out)
}
