library(cmdstanr);library(here);library(rlist)

#load in species observations as a list
fish_reef_trim_3403_2=read.csv(here('data','sp_info.csv'))

rvc_occs_3403=list()
reef_occs_3403=list()
for(q in 1:87){
  rvc_occs_3403[[q]]=read.csv(here('data','RVC','species observations',paste(sprintf("%02d",q),'_RVC_',gsub(' ', '_',fish_reef_trim_3403_2$commonname[q]),'.csv',sep='')))
  reef_occs_3403[[q]]=read.csv(here('data','REEf','species observations',paste(sprintf("%02d",q),'_REEF_',gsub(' ', '_',fish_reef_trim_3403_2$commonname[q]),'.csv',sep='')))
}
save.image(here('data','data.RData'))
####Stan models####
set_cmdstan_path()

file_mv= file.path(cmdstan_path(),"multi_comb_survey_mod.stan")
mod_mv=cmdstanr::cmdstan_model(file_mv)
mod_diag=data.frame(sp=fish_reef_trim_3403_2$commonname,rhat.flag=NA,ess.bulk.flag=NA,ess.tail.flag=NA)  

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
     iter_sampling = 800,
     refresh = 100,
     max_treedepth = 20 # print update every 500 iters
   )

   params<- fitmv$draws(format='df')
   
   plot_path<- here('outputs','figures')
   scaled_mv_timeseries_plot(q=q,ts1=rvc_ts_3403[[q]],ts2=reef_ts_3403[[q]],sp=fish_reef_trim_3403_2$commonname[q],params=params,path=plot_path,TT=26,TT.rvc=23,yr.start=1993,yr.end=2018)

   mod_par_path<- here('outputs','parameter estimates')
   pars = c('cut','a_hab1','a_hab2','a_strat1','a_strat2','sd_strat2','sd_strat1','sd_psu','sd_hab1','sd_hab2','sd_mth1','sd_mth2','sd_site','sd_dv','sd_dmy','sd_my','sd_r','sd_q','x','a_yr1','a_yr2','beta1','beta2','phi','Cor_t')
    
   write.csv(fitmv$draws(format='df',variables=pars),file.path(mod_par_path,paste(sprintf("%02d",q),'_',gsub(' ', '_',fish_reef_trim_3403_2$commonname[q]),'.csv',sep='')))
   
   mod_par_path_2<- here('outputs','species parameter summary')
   
   write.csv(as.data.frame(fitmv$summary(variables=pars)),file.path(mod_par_path_2,paste(sprintf("%02d",q),'_',gsub(' ', '_',fish_reef_trim_3403_2$commonname[q]),'_model1','.csv',sep='')))
  
   #posterior predictive check plots
   yrep1=fitmv$draws(variables='y_rep_rvc')
   yrep2=fitmv$draws(variables='y_rep_reef')
   
   ppc_comp_plot(emp1=spp_rvc,emp2=spp_reef,yrep1=yrep1,yrep2=yrep2,pdf=T,file.name = paste(sprintf("%02d",q),'_',gsub(' ', '_',fish_reef_trim_3403_2$commonname[q]),sep=''),path=here('outputs','model diagnostics','ppc'))
   
   #model diagnostic check
   pars_all = c('cut','a_hab1','a_hab2','a_strat1','a_strat2','sd_strat2','sd_strat1','sd_psu','sd_hab1','sd_hab2','sd_mth1','sd_mth2','sd_site','sd_dv','sd_dmy','sd_my','sd_r','sd_q','x','a_yr1','a_yr2','beta1','beta2','a_mth1','a_mth2','a_psu','a_dmy','a_my','phi','Cor_t')
   
   rhat=fitmv$summary(pars_all)$rhat
   mod_diag[q,2]=any(na.omit(rhat)>=1.05) #rhat = 1.05 cut-off
   ess1=fitmv$summary(pars_all)$ess_bulk
   ess2=fitmv$summary(pars_all)$ess_tail
   mod_diag[q,3]=any(na.omit(ess1)<400) #n= 400 min. effective sample size cut-off 
   mod_diag[q,4]=any(na.omit(ess2)<400)
   write.csv(mod_diag[q,],here('outputs','model diagnostics',paste(sprintf("%02d",q),'_',fish_reef_trim_3403_2$commonname[q],'_diag.csv',sep='')))
}
write.csv(mod_diag,here('outputs','model diagnostics','model_diag_summary.csv'))

###Survey effort sensitivity test####
#for 9 sp with highest agreement
sp9=c(15,22,60,37,14,66,24,11,44)

rho_sens=list()
for(z in 2:length(sp10)){
  spp_rvc<- rvc_occs_3403[[sp5[z]]]
  spp_reef<- reef_occs_3403[[sp5[z]]]
  
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
  
  #sub-sample down to different fractions by year
  #75%, 50%, 25% sample effort
  n_s_reef =spp_reef %>% group_by(year) %>% summarize(n=length(unique(formid)))
  n_s_reef$n_sub75=round(n_s_reef$n*0.75)
  n_s_reef$n_sub50=round(n_s_reef$n*0.5)
  n_s_reef$n_sub25=round(n_s_reef$n*0.25)
  
  spp_reef_25=list()
  data_25=list()
  rho_sens[[z]]=matrix(nrow=1800,ncol=20)
  for(q in 1:20){ #randomly sample observations within each year to reduce overall sampling effort by 25, 50, 75%
   # spp_reef_75[[q]]=subset(spp_reef,year==min(year))[sample(nrow(subset(spp_reef,year==min(year))),size=n_s_reef$n_sub75[1],replace=FALSE),]
#    spp_reef_50[[q]]=subset(spp_reef,year==min(year))[sample(nrow(subset(spp_reef,year==min(year))),size=n_s_reef$n_sub50[1],replace=FALSE),]
    spp_reef_25[[q]]=subset(spp_reef,year==min(year))[sample(nrow(subset(spp_reef,year==min(year))),size=n_s_reef$n_sub25[1],replace=FALSE),]
    for(t in 2:nrow(n_s_reef)){
   #   spp_reef_75[[q]]=rbind(spp_reef_75[[q]],subset(spp_reef,year==n_s_reef$year[t])[sample(n_s_reef$n[t],size=n_s_reef$n_sub75[t],replace=FALSE),])
 #     spp_reef_50[[q]]=rbind(spp_reef_50[[q]],subset(spp_reef,year==n_s_reef$year[t])[sample(n_s_reef$n[t],size=n_s_reef$n_sub50[t],replace=FALSE),])
      spp_reef_25[[q]]=rbind(spp_reef_25[[q]],subset(spp_reef,year==n_s_reef$year[t])[sample(n_s_reef$n[t],size=n_s_reef$n_sub25[t],replace=FALSE),])
    }
  
    X1<- matrix(data=c(scale(as.numeric(spp_rvc$DEPTH)),scale(as.numeric(spp_rvc$DEPTH)^2)),ncol=2,nrow=nrow(spp_rvc))
  X2<- matrix(data=c(scale(as.numeric(spp_reef_25[[q]]$btime)),scale(as.numeric(spp_reef_25[[q]]$averagedepth)),scale(as.numeric(spp_reef_25[[q]]$averagedepth)^2),scale(as.numeric(spp_reef_25[[q]]$visibility)),scale(as.numeric(spp_reef_25[[q]]$current)),spp_reef_25[[q]]$exp_binary),ncol=6,nrow=nrow(spp_reef_25[[q]]))

  data_25[[q]]=list(y1 = spp_rvc$NUM.total2,
                    y2 = spp_reef_25[[q]]$abundance2,
                    N1 = nrow(spp_rvc),
                    N2 = nrow(spp_reef_25[[q]]),
                    N_psu = length(unique(spp_rvc$PSU_YEAR)),
                    psu_yr = as.numeric(factor(spp_rvc$PSU_YEAR)),
                    N_hab1 = length(unique(spp_rvc$HABITAT_CD)),
                    hab_class1=as.numeric(factor(spp_rvc$HABITAT_CD)),
                    N_strat1=length(unique(spp_rvc$STRAT)),
                    stratum1=as.numeric(factor(spp_rvc$STRAT)),
                    N_mth1=length(unique(spp_rvc$MONTH)),
                    mth1=as.numeric(factor(spp_rvc$MONTH)),
                    N_hab2 = length(unique(spp_reef_25[[q]]$hab_class)),
                    hab_class2=as.numeric(factor(spp_reef_25[[q]]$hab_class)),
                    site=as.numeric(factor(spp_reef_25[[q]]$geogr)),
                    N_site=length(unique(spp_reef_25[[q]]$geogr)),
                    N_strat2=length(unique(spp_reef_25[[q]]$stratum)),
                    stratum2=as.numeric(factor(spp_reef_25[[q]]$stratum)),
                    N_mth2=length(unique(spp_reef_25[[q]]$month)),
                    mth2=as.numeric(factor(spp_reef_25[[q]]$month)),
                    diver=as.numeric(factor(spp_reef_25[[q]]$fish_memberid)),
                    N_dv=length(unique(spp_reef_25[[q]]$fish_memberid)),
                    dmy=as.numeric(factor(spp_reef_25[[q]]$site_dmy)),
                    N_dmy=length(unique(spp_reef_25[[q]]$site_dmy)),
                    my=as.numeric(factor(spp_reef_25[[q]]$mth_cluster)),
                    N_my=length(unique(spp_reef_25[[q]]$mth_cluster)),
                    K=max(spp_reef_25[[q]]$abundance2),
                    X1=X1,
                    Z1=ncol(X1),
                    X2=X2,
                    Z2=ncol(X2),
                    TT=max(spp_reef_25[[q]]$year)-min(spp_reef_25[[q]]$year)+1,
                    N_yr1=length(unique(spp_rvc$YEAR)),
                    yr_index1=unique((spp_rvc$YEAR-min(spp_rvc$YEAR))+1),
                    year_id1=(spp_rvc$YEAR-min(spp_rvc$YEAR))+1,
                    N_yr2=length(unique(spp_reef_25[[q]]$year)),
                    yr_index2=sort(unique(as.numeric(factor(spp_reef_25[[q]]$year)))),
                    year_id2=as.numeric(factor(spp_reef_25[[q]]$year))
  )
  
  fit_mv25=  mod_mv$sample(
    data = data_25[[q]],
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 200,
    iter_sampling = 300,
    refresh = 100,
    max_treedepth = 20 # print update every 500 iters
  )
    

  rho_sens[[z]][,q]=as.data.frame(fit_mv25$draws(variable='Cor_t[1,2]',format='draws_matrix'))[,1]
  }
  
  #output rho for each subset across the 20 random sets:
  write.csv(rho_sens[[z]],here('outputs','sampling effort sensitivity',paste(fish_reef_trim_3403_2$commonname[sp5[z]],'sens_n.csv',sep='')))
}


#Independent index fits####
#independent fits
file_sep= file.path(cmdstan_path(),"mod2.stan")
mod_sep=cmdstanr::cmdstan_model(file_sep)

#fit uncorrelated state-space models to each survey
#outputs
#figures of fits
#parameter summaries
#parameter estimates

#mv fit vs ind fit
for(q in 1:length(rvc_occs_3403)){
  spp_rvc<- rvc_occs_3403[[sp10[q]]]
  spp_reef<- reef_occs_3403[[sp10[q]]]
  
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
            yr_index1=sort(unique(spp_rvc$YEAR)-min(spp_rvc$YEAR)+1),
            year_id1=as.numeric(factor(spp_rvc$YEAR)),
            N_yr2=length(unique(spp_reef$year)),
            yr_index2=sort(unique(as.numeric(factor(spp_reef$year)))),
            year_id2=as.numeric(factor(spp_reef$year))
  )
  
  fitsep<- mod_sep$sample(
    data = data,
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 200,
    iter_sampling = 500,
    refresh = 100,
    max_treedepth = 20 # print update every 500 iters
  )
 
   params<- as.data.frame(fitsep$draws(format='df'))
  
  plot_path<- here('outputs','figures','id fits')
  scaled_timeseries_plot(i=q,ts1=rvc_ts_3403[[sp10[q]]],ts2=reef_ts_3403[[sp10[q]]],mod='model2',sp=fish_reef_trim_3403_2$commonname[sp10[q]],params1=params,params2=params,path=plot_path,TT=26,TT.rvc=23,yr.start=1993,yr.end=2018,n.iter=3000)
  
  mod_par_path<- here('outputs','parameter estimates','id fits')
  pars = c('cut','a_hab1','a_hab2','a_strat1','a_strat2','sd_strat2','sd_strat1','sd_psu','sd_hab1','sd_hab2','sd_mth1','sd_mth2','sd_site','sd_dv','sd_dmy','sd_my','sd_r1','sd_q1','sd_r2','sd_q2','x1','x2','a_yr1','a_yr2','beta1','beta2','phi')
  
  write.csv(fitsep$draws(format='df',variables=pars),file.path(mod_par_path,paste(sprintf("%02d",sp10[q]),'_',gsub(' ', '_',fish_reef_trim_3403_2$commonname[sp10[q]]),'.csv',sep='')))
  
  mod_par_path_2<- here('outputs','species parameter summary', 'id fits')
  
  write.csv(as.data.frame(fitsep$summary(variables=pars)),file.path(mod_par_path_2,paste(sprintf("%02d",sp10[q]),'_',gsub(' ', '_',fish_reef_trim_3403_2$commonname[sp10[q]]),'_model2','.csv',sep='')))

}

 