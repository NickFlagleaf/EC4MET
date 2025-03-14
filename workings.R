env.info<-read.csv("C:/Users/fra338/OneDrive - CSIRO/Documents/PDF project/CAIGE data/All_environments.csv")

library(devtools)
library(roxygen2)
library(usethis)
library(FunFunctions)

#Source functions----
source("C:/Users/fra338/OneDrive - CSIRO/Documents/PDF project/CAIGE data/GxE analysis/AS REML R analyses/2 stage analysis/FA.nxt.SVs function.R")
source("C:/Users/fra338/OneDrive - CSIRO/Documents/PDF project/CAIGE data/GxE analysis/AS REML R analyses/2 stage analysis/fa.sum function.R")

load("~/PDF project/CAIGE data/GxE analysis/AS REML R analyses/All genomic A add only fitted FA models 2011 to 2024 +newlines preds.RData")
all.mods<-list(DIAGDAIG=lmmod.diag,fa1fa1=lmmod.fa1.fa1,fa2fa1=lmmod.fa2.fa1,fa3fa1=lmmod.fa3.fa1)
lapply(all.mods,function(x) x$converge)
all.fa.sums<-lapply(all.mods[-1],function(x) fa.sum(x,addGfac = "GID",non.add = F,add = T,Efac = "Env"))


env.info<-env.info[env.info$Year%in%c(2020:2022),c("Environment","Lat","Long","Year","Sowing.date")]
env.info<-env.info[!env.info$Sowing.date=="",]
env.info<-env.info[env.info$Environment %in% names(all.fa.sums$fa3fa1$Fixed.effects),]

env.info<-cbind(env.info,
                "Main_E_effect"=all.fa.sums$fa3fa1$Fixed.effects[env.info$Environment],
                all.fa.sums$fa3fa1$Additive_G_effects$Lam[env.info$Environment,])
CAIGE20_22envs<-env.info
save(CAIGE20_22envs,file="data/CAIGE20_22envs.RData")



wether.var.cors<-sapply(1:ncol(weather.ECs.SILO),function(x) cor(weather.ECs.BARRA[,x],weather.ECs.SILO[,x]))
names(wether.var.cors)<-colnames(weather.ECs.SILO)


fullcoords<-expand.grid("lon"=seq(110,160,length.out=200),"lat"=seq(-20,-40,length.out=120))
wheatarearast<-raster::raster(x = "C:/Users/fra338/OneDrive - CSIRO/Documents/PDF project/Climate projections/CMIP6/monfreda/wheat_HarvestedAreaFraction.tif")
wheat.areaspergrid<-terra::extract(wheatarearast,y = fullcoords)
wheat.areaspergrid<-cbind(fullcoords,wheat.areaspergrid)
wheat.areaspergrid<-wheat.areaspergrid[wheat.areaspergrid$wheat.areaspergrid>0.02 & !is.na(wheat.areaspergrid$wheat.areaspergrid),]

locs<-paste(round(wheat.areaspergrid$lon,3),round(wheat.areaspergrid$lat,3),sep="_")
locs<-paste(rep(2018:2020,each=nrow(wheat.areaspergrid)),rep(locs,3),sep = "_")
wheat.area.envs<-data.frame("Env"=locs,
                            "Year"=sapply(locs,function(x) strsplit(x,"_")[[1]][1]),
                            "Lon"=sapply(locs,function(x) strsplit(x,"_")[[1]][2]),
                            "Lat"=sapply(locs,function(x) strsplit(x,"_")[[1]][3]))
wheat.area.envs$sow.dates<-paste("01/05/",wheat.area.envs$Year,sep="")

save(wheat.area.envs,file="data/wheat.area.envs.RData")

2018_148.442_-20.168,2018_151.457_-28.908,2018_115.779_-29.58,2018_150.955_-29.916 2018_117.286_-33.445 2018_134.874_-33.445 2018_147.437_-33.613 2018_140.402_-34.286 2018_139.648_-34.622 2018_136.131_-34.958 2018_146.181_-35.294 2019_148.442_-20.168 2019_151.457_-28.908 2019_115.779_-29.58 2019_150.955_-29.916 2019_117.286_-33.445 2019_134.874_-33.445 2019_147.437_-33.613 2019_140.402_-34.286 2019_139.648_-34.622 2019_136.131_-34.958 2019_146.181_-35.294 2020_148.442_-20.168 2020_151.457_-28.908 2020_115.779_-29.58 2020_150.955_-29.916 2020_117.286_-33.445 2020_134.874_-33.445 2020_147.437_-33.613 2020_140.402_-34.286 2020_139.648_-34.622 2020_136.131_-34.958 2020_146.181_-35.294

data("CAIGE20_22envs")
SILO.wthr<-get.SILO.weather(Envs = CAIGE22_23envs$Environment,
                            Lats = CAIGE22_23envs$Lat,
                            Lons = CAIGE22_23envs$Long,
                            Years = CAIGE22_23envs$Year)

BARRA.wthr<-get.BARRA.weather(Envs = CAIGE22_23envs$Environment,
                              Lats = CAIGE22_23envs$Lat,
                              Lons = CAIGE22_23envs$Long,
                              Years = CAIGE22_23envs$Year)

data("wheat.area.envs")
wheat.area.wthr<-get.SILO.weather(Envs = wheat.area.envs$Env,
                                  Lats = wheat.area.envs$Lat,
                                  Lons = wheat.area.envs$Lon,
                                  Years = wheat.area.envs$Year)


#Get weather ECs
weather.ECs.SILO<-get.W.ECs(weather = SILO.wthr,
                            sow.dates = CAIGE22_23envs$Sowing.date)

weather.ECs.BARRA<-get.W.ECs(weather = BARRA.wthr,
                            sow.dates = CAIGE22_23envs$Sowing.date)

wheat.area.weather.ECs.SILO<-get.W.ECs(weather = wheat.area.wthr,
                                       sow.dates = wheat.area.envs$sow.dates)

#Get soil ECs
soil.ECs<-get.S.ECs(Envs = CAIGE20_22envs$Environment,
                    Lats = CAIGE20_22envs$Lat,
                    Lons = CAIGE20_22envs$Long)

wheat.area.soil.ECs<-get.S.ECs(Envs = wheat.area.envs$Env,
                               Lats = wheat.area.envs$Lat,
                               Lons = wheat.area.envs$Lon)

#Combine weather and soil data
obsECs<-cbind(weather.ECs.SILO$,soil.ECs)
wheat.area.ECs<-cbind(wheat.area.weather.ECs.SILO,wheat.area.soil.ECs)


obs.env.effs<-CAIGE22_23envs[,c("Main_E_effect","FA1","FA2","FA3")]
rownames(obs.env.effs)<-CAIGE22_23envs$Environment
                             
wheat.area.preds<-pred.env.effs(train.ECs=obsECs,
                                new.ECs=wheat.area.ECs,
                                E.effs=obs.env.effs)
                      

 
#Work out iclasses
iclasses<-apply(wheat.area.preds[,c("FA1","FA2","FA3")],1,function(x) paste(ifelse(x>0,"p","n"),collapse = ""))

#plot map    
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

{
stime<-Sys.time()
wheat.area.wthr<-get.SILO.weather(Envs = wheat.area.envs$Env[wheat.area.envs$Year=="2018"],
                                  Lats = wheat.area.envs$Lat[wheat.area.envs$Year=="2018"],
                                  Lons = wheat.area.envs$Lon[wheat.area.envs$Year=="2018"],
                                  Years = wheat.area.envs$Year[wheat.area.envs$Year=="2018"])
end.time<-Sys.time()
end.time-stime
}








#Test CMPI6 function
cmip.weather<-get.CMIP6.weather(Envs = CAIGE20_22envs$Environment,
                                Lats = CAIGE20_22envs$Lat,
                                Lons = CAIGE20_22envs$Long,
                                Years = CAIGE20_22envs$Year+18,
                                SSPs = c("ssp585","ssp370","ssp245","ssp126"),
                                GCMs = c("ACCESS-CM2","ACCESS-ESM1-5","CESM2","CMCC-ESM2","CNRM-ESM2-1",
                                         "EC-Earth3","MPI-ESM1-2-HR","NorESM2-MM","UKESM1-0-LL"),
                                ncores = 3)





