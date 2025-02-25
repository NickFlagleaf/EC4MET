


#Load example dataset of CAIGE 2023 and 2024 environments
data("CAIGE23_24envs")

#Get daily weather data from SILO for each environment
wthr<-Get.SILO.weather(Envs = CAIGE23_24envs$Environment,
                       Lats = CAIGE23_24envs$Lat,
                       Lons = CAIGE23_24envs$Long,
                       Years = CAIGE23_24envs$Year)

#Calculate ECs based on the daily weather data
weatherECs<-Get.ECs(weather = wthr)
