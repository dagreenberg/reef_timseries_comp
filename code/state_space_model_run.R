#required packages
library(cmdstanr);library(loo);library(bayesplot);library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
library(here)
here()

###Functions
source(here('code','functions.R'))

####Data Set-up####
###REEF survey data
R<- read.csv(here('data','REEF','REEF_keys_data.csv')) #Reef sightings 

###REEF geographic site data
library(sp)
reef_geog<- read.csv(here('data','REEF','TWAgeog.csv'),na.strings=c('NULL'))
site_surveys<- R %>% group_by(geogr) %>% summarize(n.surv=n_distinct(formid))
reef_geog$no.surveys<- site_surveys$n.surv[match(reef_geog$geogid,site_surveys$geogr)]
reef_geog$lat_full<- reef_geog$lat #Copy of the full latitude - these are in degrees, minutes (with seconds as a fraction of minutes)
reef_geog$lon_full<- reef_geog$lon #Copy of the full longitude

reef_geog<-reef_geog%>% #Separate out degrees and minutes
  tidyr::separate(lat,into=c("lat_deg","lat_min"),sep=" ")%>%
  tidyr::separate(lon,into=c("lon_deg","lon_min"),sep=" ")

#converting to numeric
col.num<-c("lat_deg","lat_min","lon_deg","lon_min")
reef_geog[col.num]<-sapply(reef_geog[col.num],as.numeric) #will get warnings for aberrant entries

#converting from degrees minutes seconds to decimal degrees
reef_geog$lat_dd<- reef_geog$lat_deg+reef_geog$lat_min/60
reef_geog$lon_dd<- reef_geog$lon_deg-reef_geog$lon_min/60
reef_geog<- reef_geog[complete.cases(reef_geog$lat_dd),] #remove sites/regions without spatial coordinates
reef_geog<- reef_geog[complete.cases(reef_geog$lon_dd),] #remove sites/regions without spatial coordinates

R$site4=reef_geog$region.id[match(R$site,reef_geog$geogid)]
###Reef Visual Census (RVC) survey data
fk_99_18<- read.csv(here('data','RVC','Florida_Keys_RVC.csv')) #RVC survey data from 1999 to 2018
fk_79_98<- read.csv(here('data','RVC','Pre1999_edit.txt')) #RVC survey dating from 1979 to 1998

#Synonymize column names between the two datasets
colnames(fk_79_98)<- c(colnames(fk_99_18)[1],colnames(fk_99_18)[3:5],colnames(fk_99_18)[2],colnames(fk_99_18)[6],colnames(fk_99_18)[7:18],colnames(fk_99_18)[19],colnames(fk_99_18)[20]) 
#Join together datasets
fk_79_18<- full_join(fk_99_18,fk_79_98)

#Identify distinct RVC sampling sites
fk_79_18$LAT_LON<- paste(fk_79_18$LAT_DEGREES,fk_79_18$LON_DEGREES,sep='_') #Find unique geographic position
rvc_sites<- data.frame(lat=fk_79_18$LAT_DEGREES,lon=fk_79_18$LON_DEGREES,lat_lon=fk_79_18$LAT_LON) 
rvc_sites<- distinct(rvc_sites,lat_lon,.keep_all = T) #keep unique sample sites from the RVC surveys

#Identify the closest REEF dive site for each RVC sampling unit (PSU)
library(rgeos)
set1 <- SpatialPoints(cbind(rvc_sites$lat,rvc_sites$lon)) #Set of RVC site points
set2 <- SpatialPoints(cbind(reef_geog$lat_dd,reef_geog$lon_dd)) #Set of REEF site points
#set3 <- SpatialPoints(cbind(rvc_grid$lat,rvc_grid$lon)) #Set of REEF site points
matched_reef_site<- apply(gDistance(set1, set2, byid=TRUE), 2, which.min)
#matched_grid<- apply(gDistance(set1, set3, byid=TRUE), 2, which.min)#Find closest REEF site by geographic distance between point sets
rvc_sites$REEF_site <- reef_geog$geog[matched_reef_site]
rvc_sites$REEF_lat<- reef_geog$lat_dd[matched_reef_site]
rvc_sites$REEF_lon<- reef_geog$lon_dd[matched_reef_site]
rvc_sites$region.id<- reef_geog$region.id[matched_reef_site]
fk_79_18[,25:28]<- rvc_sites[match(fk_79_18$LAT_LON,rvc_sites$lat_lon),4:7]

###Geographic filtering 
#Intersect REEF with the RVC benthic habitat map
library(sf)
rvc_grid<-st_read(dsn=here('data','RVC Grid'), layer='FlaKeys_Grid') #Read in Florida Keys benthic habitat sampling grid
#rvc_grid<- st_transform(rvc_grid,4326) #Set to WGS84

#3403 regions
reef_pts<-st_as_sf(x = reef_geog, 
                   coords = c("lon_dd", "lat_dd"),
                   crs = 4326)
reef_pts<- st_transform(reef_pts, "+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83")
rvc_pts<-st_as_sf(x = rvc_sites, 
                  coords = c("lon", "lat"),
                  crs = 4326)
rvc_pts<- st_transform(rvc_pts, "+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83")

grid_match<- st_intersects(reef_pts,rvc_grid,sparse=T)
reef_geog$grid_match<- NA
for(i in 1:nrow(reef_pts)){
  if(length(grid_match[[i]])==1){
    reef_geog$grid_match[i]=grid_match[[i]]  
  }else{
    reef_geog$grid_match[i]=NA
  }
}
reef_geog$hab_class<- as.factor(rvc_grid$habclass[reef_geog$grid_match]) #Extracts the habitat class from the sampling grid for each REEF dive site
reef_geog$stratum<- as.factor(rvc_grid$STRAT[reef_geog$grid_match]) #Extracts the habitat class from the sampling grid for each REEF dive site
reef_matched<- subset(reef_geog,is.na(grid_match)==F)

grid_match_rvc<- st_intersects(rvc_pts,rvc_grid,sparse=T) #This matches the RVC survey data to the sampling grid to make sure benthic habitat classes are accurate for the older survey data

#test=apply(st_distance(rvc_pts,rvc_grid, by_element = T, which = "Euclidean"), 2, which.min)
#closest_grid<- apply(gDistance(rvc_pts,rvc_grid, byid=TRUE), 2, which.min)
rvc_sites$grid_match<- NA
for(i in 1:nrow(rvc_sites)){
  if(length(grid_match_rvc[[i]])==1){
    rvc_sites$grid_match[i]=grid_match_rvc[[i]]  
  }else{
    rvc_sites$grid_match[i]=NA
  }
}
rvc_sites$hab_class<- as.factor(rvc_grid$habclass[rvc_sites$grid_match])
rvc_sites$stratum<- as.factor(rvc_grid$STRAT[rvc_sites$grid_match])
rvc_sites$utm_x<- rvc_grid$utm_x[rvc_sites$grid_match]
rvc_sites$utm_y<- rvc_grid$utm_y[rvc_sites$grid_match]
rvc_sites1<- subset(rvc_sites, is.na(grid_match)==F)#remove sampling plots outside of the grid & removed plots in inshore patch reefs (low number of samples and not sampled by REEF)
rvc_sites2<- subset(rvc_sites, is.na(grid_match)==T)#remove sampling plots outside of the grid & removed plots in inshore patch reefs (low number of samples and not sampled by REEF)
closest_grid<- apply(gDistance(SpatialPoints(cbind(rvc_sites2$lat,rvc_sites2$lon)),SpatialPoints(cbind(rvc_grid$lat,rvc_grid$lon)), byid=TRUE), 2, which.min)
rvc_sites2$hab_class<- rvc_grid$habclass[closest_grid]
rvc_sites2$stratum<- rvc_grid$STRAT[closest_grid]
rvc_sites<- rbind(rvc_sites1,rvc_sites2)

###Data filtering and processing (see functions)
#Fish data from REEF - remove ultra rare and basket species designations
fish_reef<- read.csv(here('data','REEF',"Caribbean_fish_trait_matrix.csv")) #fish species found in the Tropical Western Atlantic
fish_rvc<- read.csv(here('data','RVC',"Florida_keys_taxonomic_data.csv"))
fish_rvc<- subset(fish_rvc,gsub('.*\\ ', '', fish_rvc$SCINAME)!='sp.') #remove unknown species
fish_rvc<- subset(fish_rvc, SCINAME %in% fish_reef$sciname2)
m<- match(fish_reef$sciname2,fish_rvc$SCINAME)
fish_reef$rvc_code<- fish_rvc$SPECIES_CD[m]

fk_93_18<- subset(fk_79_18,YEAR>=1993) #Subset RVCdataset from 1993 to match the first year of REEF surveys

###Key Largo data set-up####
rvc_occs_3403_1<- rvc_filter(fk_93_18,GZ='3403',sites=rvc_sites,sp=fish_reef,year.start = 1993,year.end = 2018,strat.omit = 'n',hab.omit = 'n')
rvc_ts_3403_1<- ts_rvc(rvc_occs_3403_1,miss="T")
rvc_ts_filter_3403<- rlist::list.filter(rvc_ts_3403_1,length(na.omit(p.occ))>18) #Removes all species that had >30% of years where it was not sighted in the RVC dataset

for(i in 1:nrow(fish_reef)){
  fish_reef$rvc_sight_freq_3403[i]<- mean(na.omit(rvc_ts_3403_1[[i]]$p.occ))
}
rvc_green_3403<- do.call(rbind, lapply(rvc_ts_filter_3403, data.frame, stringsAsFactors=FALSE))
rvc_green_sp_3403<- unique(rvc_green_3403$sp)
fish_reef_trim_3403<- subset(fish_reef, rvc_code %in% rvc_green_sp_3403)


reef_occs_3403_1<- reef_filter(R,GZ='3403',sp=fish_reef_trim_3403,geog=reef_matched,year.start = 1993,year.end = 2018,strat.omit = 'n',hab.omit = 'n')
reef_ts_3403_1<- ts_reef(reef_occs_3403_1,sp=fish_reef_trim_3403)

for(i in 1:nrow(fish_reef_trim_3403)){
  fish_reef_trim_3403$reef_sight_freq_3403[i]<- mean(na.omit(reef_ts_3403_1[[i]]$p.occ))
}

fish_reef_trim_3403<- subset(fish_reef_trim_3403, rvc_sight_freq_3403>0.01) #Only keep those species with sufficient data

##Species to drop - reef vagrants, cave dwellers, taxonomically challenging or recently identified splits 
drop<- c("Scrawled Cowfish","Purple Reeffish","Beaugregory","Dusky Damselfish",
         "Longfin Damselfish","Reef Croaker","Black Margate",
         "White Margate","Blue Runner","Yellow Jack","Cero","Bluelip Parrotfish",
         "Bucktooth Parrotfish","Green Razorfish","Rosy Razorfish","Tobaccofish","Sharksucker",
         "Glassy Sweeper","Atlantic Spadefish","Yellowhead Jawfish","Cubbyu","Highhat",
         "Tarpon","Golden Smooth Trunkfish") #Most of these are too rare in the surveys but not initially filtered out

fish_reef_trim_3403_2<- subset(fish_reef_trim_3403,commonname %notin% drop)
rownames(fish_reef_trim_3403_2)<- seq(1:nrow(fish_reef_trim_3403_2))
fish_reef_trim_3403_2<- fish_reef_trim_3403_2[complete.cases(fish_reef_trim_3403_2),]
write.csv(fish_reef_trim_3403_2,here('outputs','sp_info.csv'))

rvc_occs_3403<- rvc_filter(fk_93_18,GZ='3403',sites=rvc_sites,sp=fish_reef_trim_3403_2,year.start = 1993,year.end = 2018,strat.omit = 'n',hab.omit = 'n')
reef_occs_3403<- reef_filter(R,GZ='3403',sp=fish_reef_trim_3403_2,geog=reef_matched,year.start = 1993,year.end = 2018,strat.omit = 'n',hab.omit = 'n')
reef_ts_3403<- ts_reef(reef_occs_3403,sp=fish_reef_trim_3403_2)
rvc_ts_3403<- ts_rvc(rvc_occs_3403,miss='F') #Turns

#Final sites included in filtered datasets:
rvc_geog_3403<- rvc_occs_3403[[1]] %>% distinct(LAT_LON,.keep_all = T)
write.csv(rvc_geog_3403,file.path(here('outputs','geographic data'),'rvc_key_largo_sites.csv'))
reef_geog_3403<- reef_occs_3403[[1]] %>% group_by(site) %>% summarize(n=n_distinct(formid),hab=unique(hab_class2),strat=unique(stratum),lat.dd=unique(lat_dd),lon.dd=unique(lon_dd))
write.csv(reef_geog_3403,file.path(here('outputs','geographic data'),'reef_key_largo_sites.csv'))

save.image(here('data','data.RData'))


####Stan models####
set_cmdstan_path(here('code','stan models'))

file_mv= file.path(cmdstan_path(),"multi_comb_survey_mod.stan")
mod_mv=cmdstanr::cmdstan_model(file_mv)
  
#Multivariate model fit to every species####
 for(q in 1:length(rvc_occs_3403)){
   spp_rvc<- rvc_occs_3403[[q]]
   spp_reef<- reef_occs_3403[[q]]
   
   #Turn model factors into numeric indices for Stan
   spp_rvc$psu_yr = as.numeric(factor(spp_rvc$PSU_YEAR))
   spp_rvc$hab1=as.numeric(factor(spp_rvc$HAB_CD2))
   spp_rvc$stratum1=as.numeric(factor(spp_rvc$STRAT))
   spp_rvc$mth1=as.numeric(factor(spp_rvc$MONTH))
   
   spp_reef$hab2=as.numeric(factor(spp_reef$hab_class2))
   spp_reef$site=as.numeric(factor(spp_reef$geogr))
   spp_reef$stratum2=as.numeric(factor(spp_reef$stratum))
   spp_reef$mth2=as.numeric(factor(spp_reef$month))
   spp_reef$diver=as.numeric(factor(spp_reef$fish_memberid))
   spp_reef$dmy=dmy=as.numeric(factor(spp_reef$site_dmy))
   spp_reef$my=dmy=as.numeric(factor(spp_reef$mth_cluster))
   
   #matrix for covariates estimated with fixed effects
   X1<- matrix(data=c(scale(as.numeric(spp_rvc$DEPTH)),scale(as.numeric(spp_rvc$DEPTH)^2)),ncol=2,nrow=nrow(spp_rvc))
   X2<- matrix(data=c(scale(as.numeric(spp_reef$btime)),scale(as.numeric(spp_reef$averagedepth)),scale(as.numeric(spp_reef$averagedepth)^2),scale(as.numeric(spp_reef$visibility)),scale(as.numeric(spp_reef$current)),spp_reef$exp_binary),ncol=6,nrow=nrow(spp_reef))

   data=list(y1 = spp_rvc$NUM.total2,
             y2 = spp_reef$abundance2,
             N1 = nrow(spp_rvc),
             N2 = nrow(spp_reef),
             N_psu = length(unique(spp_rvc$PSU_YEAR)),
             psu_yr = as.numeric(factor(spp_rvc$PSU_YEAR)),
             N_hab1 = length(unique(spp_rvc$HABITAT_CD)),
             hab_class1=as.numeric(factor(spp_rvc$HABITAT_CD)),
             N_strat1=length(unique(spp_rvc$STRAT)),
             stratum1=as.numeric(factor(spp_rvc$STRAT)),
             N_mth1=length(unique(spp_rvc$MONTH)),
             mth1=as.numeric(factor(spp_rvc$MONTH)),
             N_hab2 = length(unique(spp_reef$hab_class)),
             hab_class2=as.numeric(factor(spp_reef$hab_class)),
             site=as.numeric(factor(spp_reef$geogr)),
             N_site=length(unique(spp_reef$geogr)),
             N_strat2=length(unique(spp_reef$stratum)),
             stratum2=as.numeric(factor(spp_reef$stratum)),
             N_mth2=length(unique(spp_reef$month)),
             mth2=as.numeric(factor(spp_reef$month)),
             diver=as.numeric(factor(spp_reef$fish_memberid)),
             N_dv=length(unique(spp_reef$fish_memberid)),
             dmy=as.numeric(factor(spp_reef$site_dmy)),
             N_dmy=length(unique(spp_reef$site_dmy)),
             my=as.numeric(factor(spp_reef$mth_cluster)),
             N_my=length(unique(spp_reef$mth_cluster)),
             K=max(spp_reef$abundance2),
             X1=X1,
             Z1=ncol(X1),
             X2=X2,
             Z2=ncol(X2),
             TT=max(spp_reef$year)-min(spp_reef$year)+1,
             N_yr1=length(unique(spp_rvc$YEAR)),
             yr_index1=unique((spp_rvc$YEAR-min(spp_rvc$YEAR))+1),
             year_id1=(spp_rvc$YEAR-min(spp_rvc$YEAR))+1,
             N_yr2=length(unique(spp_reef$year)),
             yr_index2=sort(unique(as.numeric(factor(spp_reef$year)))),
             year_id2=as.numeric(factor(spp_reef$year))
)

   fitmv<- mod_mv$sample(
     data = data,
     chains = 6, 
     parallel_chains = 6,
     iter_warmup = 200,
     iter_sampling = 500,
     refresh = 100,
     max_treedepth = 20 # print update every 500 iters
   )

   params<- fitmv2$draws(format='df')
   
   plot_path<- here('outputs','figures')
   scaled_mv_timeseries_plot(q=q,ts1=rvc_ts_3403[[q]],ts2=reef_ts_3403[[q]],sp=fish_reef_trim_3403_2$commonname[q],params=params,path=plot_path,TT=26,TT.rvc=23,yr.start=1993,yr.end=2018)

   mod_par_path<- here('outputs','parameter estimates')
   pars = c('cut','a_hab1','a_hab2','a_strat1','a_strat2','sd_strat2','sd_strat1','sd_psu','sd_hab1','sd_hab2','sd_mth1','sd_mth2','sd_site','sd_dv','sd_dmy','sd_my','sd_r','sd_q','x','a_yr1','a_yr2','beta1','beta2','phi','Cor_t')
    
   write.csv(fitmv2$draws(format='df',variables=pars),file.path(mod_par_path,paste(sprintf("%02d",q),'_',gsub(' ', '_',fish_reef_trim_3403_2$commonname[q]),'.csv',sep='')))
   
   mod_par_path_2<- here('outputs','species parameter summary')
   
   write.csv(as.data.frame(fitmv2$summary(variables=pars)),file.path(mod_par_path_2,paste(sprintf("%02d",q),'_',gsub(' ', '_',fish_reef_trim_3403_2$commonname[q]),'_model1','.csv',sep='')))
  }


 
 