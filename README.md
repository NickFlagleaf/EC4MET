# EC4MET
## Derive ECs for METs 
Derive environmental covariates (ECs) for crop trial environments to model environmental effects and genotype by environment interactions

### Installation
To install from Github using devtools:

  ```
devtools::install_github("NickFlagleaf/EC4MET")
```

### Example work flow
Load package and example dataset of info for trial environments
  ```
library(EC4MET)
data("CAIGE22_23envs")
```
The example data set includes environment names, lat and lon values and sowing dates for 18 trial environments as well as fitted environmental main effects and 
factor loadings from a factor analytic mixed model multi-environment trial analysis

```
head(CAIGE22_23envs)
```

Get daily weather data from SILO for each environment
```
wthr<-Get.SILO.weather(Envs = CAIGE22_23envs$Environment,
                       Lats = CAIGE22_23envs$Lat,
                       Lons = CAIGE22_23envs$Long,
                       Years = CAIGE22_23envs$Year)
```

Calculate ECs based on the daily weather data
```
weather.ECs<-Get.W.ECs(weather = wthr)
```

Derive ECs based on SLGA soil data
```
soil.ECs<-Get.S.ECs(Envs = CAIGE22_23envs$Environment,
                    Lats = CAIGE22_23envs$Lat,
                    Lons = CAIGE22_23envs$Long)
```

