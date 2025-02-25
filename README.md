

# EC4MET
## Multi-trait and multi-model ensemble Genomic Selection

### Installation
To install from Github using devtools:

  ```
devtools::install_github("NickFlagleaf/EC4MET")
```


### Example work flow
Load package and example dataset of info for trial environments
  ```
library(EC4MET)
data("CAIGE23_24envs")
```
The example data set includes environment names, lat and lon values and sowing dates for 18 trial environments
```
head(CAIGE23_24envs)
```


Load example dataset of CAIGE 2023 and 2024 environments
```
data("CAIGE23_24envs")
```

Get daily weather data from SILO for each environment
```
wthr<-Get.SILO.weather(Envs = CAIGE23_24envs$Environment,
                       Lats = CAIGE23_24envs$Lat,
                       Lons = CAIGE23_24envs$Long,
                       Years = CAIGE23_24envs$Year)
```

Calculate ECs based on the daily weather data
```
weatherECs<-Get.ECs(weather = wthr)
```
