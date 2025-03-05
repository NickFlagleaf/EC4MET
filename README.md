---
output:
  html_document: default
  pdf_document: default
---
# EC4MET

***

## Derive ECs for METs 
Derive environmental covariates (ECs) for crop trial environments to model environmental effects and genotype by environment interactions

***

### Installation
To install from Github using devtools:

```
devtools::install_github("NickFlagleaf/EC4MET")
```

***

### Example work flow
Load package and example dataset of info for trial environments between 2020 and 2022:
  ```
library(EC4MET)
data("CAIGE20_22envs")
```

The example data set includes environment names, lat and lon values and sowing dates for 18 trial environments as well as fitted environmental main effects and 
factor loadings from a factor analytic mixed model multi-environment trial analysis:

```
head(CAIGE20_22envs)
```

Get daily weather data from [SILO](https://www.longpaddock.qld.gov.au/silo/) for each observed environment with the `get.SILO.weather` function:
```
obs.wthr <- get.SILO.weather(Envs = CAIGE20_22envs$Environment,
                       Lats = CAIGE20_22envs$Lat,
                       Lons = CAIGE20_22envs$Long,
                       Years = CAIGE20_22envs$Year)
```

Calculate ECs based on the daily weather data with the `get.W.ECs` function:
```
obs.weather.ECs <- get.W.ECs(weather = obs.wthr,
                           sow.dates = CAIGE20_22envs$Sowing.date)
```

Derive ECs from [SLGA](https://www.clw.csiro.au/aclep/soilandlandscapegrid/GetData-R_package.html) soil data with the `get.S.ECs` function:
```
obs.soil.ECs <- get.S.ECs(Envs = CAIGE20_22envs$Environment,
                          Lats = CAIGE20_22envs$Lat,
                          Lons = CAIGE20_22envs$Long)
```

Load example dataset for environments in the Australian grain belt from 2018 to 2020 and get weather and soil ECs for all environments.
```
data("wheat.area.envs")

wheat.area.wthr <- get.SILO.weather(Envs = wheat.area.envs$Env,
                                  Lats = wheat.area.envs$Lat,
                                  Lons = wheat.area.envs$Lon,
                                  Years = wheat.area.envs$Year)
                                  
wheat.area.wthr.ECs <- get.W.ECs(weather = wheat.area.wthr,
                               sow.dates = wheat.area.envs$sow.dates)

wheat.area.soil.ECs <- get.S.ECs(Envs = wheat.area.envs$Env,
                               Lats = wheat.area.envs$Lat,
                               Lons = wheat.area.envs$Lon)

```

Combine weather and soil ECs for the observed and all wheat belt environments 
```
obsECs <- cbind(obs.weather.ECs, obs.soil.ECs)

wheat.area.ECs <- cbind(wheat.area.wthr.ECs, wheat.area.soil.ECs)
```


Predict environmental effects for all wheat belt environments based on the observed ECs
```
obs.env.effs <- CAIGE22_23envs[ ,c("Main_E_effect","FA1","FA2","FA3")]

rownames(obs.env.effs) <- CAIGE22_23envs$Environment

wheat.area.preds <- pred.env.effs(train.ECs = obsECs,
                                new.ECs = wheat.area.ECs,
                                E.effs = obs.env.effs)
```
                                
Plot the iClasses environment types for all environments in the Australian grain belt on maps:
```
iclasses <- apply(wheat.area.preds[,c("FA1","FA2","FA3")],1,function(x) paste(ifelse(x > 0, "p", "n"), collapse = ""))

iclass.cols <- viridis::turbo(length(unique(iclasses)))
names(iclass.cols) <- unique(iclasses)
par(mfrow=c(2,2), mar=c(1,1,1,1))
  for(y in unique(wheat.area.envs$Year)){
    oz::oz()
    points(wheat.area.envs$Lon[wheat.area.envs$Year == y],
           wheat.area.envs$Lat[wheat.area.envs$Year == y], pch=15, cex=.5,
           col=iclass.cols[iclasses[wheat.area.envs$Year == y]])
    mtext(text = y, side = 3, line = -1.5, cex=1.2)
    }
```

***