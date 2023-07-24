#Output analysis
library(here); library(dplyr);library(rlist)

###Functions###
source(here('code','functions.R'))

####data load
sp_info=read.csv(here('outputs','sp_info.csv'))
dat2=read.csv(here('data','REEF','Caribbean_fish_trait_matrix.csv'))
sp_info$target=dat2$target[match(sp_info$commonname,dat2$commonname)]
sp_info$Family=dat2$Family[match(sp_info$commonname,dat2$commonname)]
pars<- list.files(here('outputs','parameter estimates'))

m_abund_list<- list() #mean abundance
r_list<- list() #temp. correlation
g_rvc_list<- list();t_rvc_kist=list()
g_reef_list<- list();t_reef_kist=list()

#Overview
sp_info$m_abund_rvc=NA;sp_info$m_abund_rvc_l90=NA;sp_info$m_abund_rvc_u90=NA;sp_info$m_abund_rvc_l95=NA;sp_info$m_abund_rvc_u95=NA
sp_info$m_abund_reef=NA;sp_info$m_abund_reef_l90=NA;sp_info$m_abund_reef_u90=NA;sp_info$m_abund_reef_l95=NA;sp_info$m_abund_reef_u95=NA
sp_info$r=NA;sp_info$r.prob=NA;sp_info$r_l80=NA;sp_info$r_u80=NA;sp_info$r_l90=NA;sp_info$r_u90=NA;sp_info$r_l95=NA;sp_info$r_u95=NA
sp_info$trend_rvc=NA;sp_info$trend_rvc_l80=NA;sp_info$trend_rvc_u80=NA;sp_info$trend_rvc_l90=NA;sp_info$trend_rvc_u90=NA;sp_info$trend_rvc_l95=NA;sp_info$trend_rvc_u95=NA
sp_info$trend_reef=NA;sp_info$trend_reef_l80=NA;sp_info$trend_reef_u80=NA;sp_info$trend_reef_l90=NA;sp_info$trend_reef_u90=NA;sp_info$trend_reef_l95=NA;sp_info$trend_reef_u95=NA
sp_info$a_hab_1.1=NA;sp_info$a_hab_1.2=NA;sp_info$a_hab_1.3=NA;sp_info$a_hab_1.4=NA;sp_info$a_hab_1.5=NA;sp_info$a_hab_1.6=NA;sp_info$a_hab_1.7=NA;sp_info$a_hab_1.8=NA;sp_info$a_hab_1.9=NA
sp_info$a_hab_2.1=NA;sp_info$a_hab_2.2=NA;sp_info$a_hab_2.3=NA;sp_info$a_hab_2.4=NA;sp_info$a_hab_2.5=NA;
sp_info$a_strat_1.1=NA;sp_info$a_strat_1.2=NA;sp_info$a_strat_1.3=NA;sp_info$a_strat_1.4=NA;sp_info$a_strat_1.5=NA;sp_info$a_strat_1.6=NA
sp_info$a_strat_2.1=NA;sp_info$a_strat_2.2=NA;sp_info$a_strat_2.3=NA;sp_info$a_strat_2.4=NA;sp_info$a_strat_2.5=NA;sp_info$a_strat_2.6=NA

for(i in 1:nrow(sp_info)){
  sp<- read.csv(here('outputs','parameter estimates',pars[i]))
  x_mat=sp[,grepl('x',colnames(sp))]
  x1=x_mat[,1:26]
  x2=x_mat[,27:ncol(x_mat)]
  c=sp[,grepl('cut',colnames(sp))]
 #mean abundance among all years:
  m_abund_rvc<- exp(apply(x1,1,mean))
  m_abund_reef<- exp(apply(log(mean_ord_to_n(x=x2,c=c)),1,mean))
  m_abund_list[[i]]<- data.frame(mean.rvc.kl=sort(m_abund_rvc),mean.reef.kl=sort(m_abund_reef),sp=sp_info$commonname[i])
  #population trend: 
  g_rvc_list[[i]]<- matrix(data=NA,nrow=3000,ncol=ncol(x1)-1)
  t_rvc_kist[[i]]<- rep(0,3000)
  for(t in 1:nrow(x1)){
    for(z in 1:(ncol(x1)-1)){
      g_rvc_list[[i]][t,z]=x1[t,z+1]-x1[t,z]
    }
    t_rvc_kist[[i]][t]<- (log(gm_mean(exp(g_rvc_list[[i]][t,]))))
  }
  
  x_abund_reef<- mean_ord_to_n(x=x2,c=c)
  g_reef_list[[i]]<- matrix(data=NA,nrow=3000,ncol=ncol(x2)-1)
  t_reef_kist[[i]]<- rep(0,3000)
  for(t in 1:nrow(x2)){
    for(z in 1:(ncol(x2)-1)){
      g_reef_list[[i]][t,z]=log(x_abund_reef[t,z+1])-log(x_abund_reef[t,z])
    }
    t_reef_kist[[i]][t]<- log(gm_mean(exp(g_reef_list[[i]][t,])))
  }
  
  r_list[[i]]=sp[,grepl('Cor_t.2.1.',colnames(sp))]
  sp_info[i,27]=median(m_abund_rvc)
  sp_info[i,28]=quantile(m_abund_rvc,0.05)
  sp_info[i,29]=quantile(m_abund_rvc,0.95)
  sp_info[i,30]=quantile(m_abund_rvc,0.025)
  sp_info[i,31]=quantile(m_abund_rvc,0.975)
  sp_info[i,32]=median(m_abund_reef)
  sp_info[i,33]=quantile(m_abund_reef,0.05)
  sp_info[i,34]=quantile(m_abund_reef,0.95)
  sp_info[i,35]=quantile(m_abund_reef,0.025)
  sp_info[i,36]=quantile(m_abund_reef,0.975)
  sp_info[i,37]=median(r_list[[i]])
  sp_info[i,38]=length(r_list[[i]][r_list[[i]]>0])/3000
  sp_info[i,39]=quantile(r_list[[i]],0.1)
  sp_info[i,40]=quantile(r_list[[i]],0.9)
  sp_info[i,41]=quantile(r_list[[i]],0.05)
  sp_info[i,42]=quantile(r_list[[i]],0.95)
  sp_info[i,43]=quantile(r_list[[i]],0.025)
  sp_info[i,44]=quantile(r_list[[i]],0.975)
  sp_info[i,45]=median(t_rvc_kist[[i]])
  sp_info[i,46]=quantile(t_rvc_kist[[i]],0.1)
  sp_info[i,47]=quantile(t_rvc_kist[[i]],0.9)
  sp_info[i,48]=quantile(t_rvc_kist[[i]],0.05)
  sp_info[i,49]=quantile(t_rvc_kist[[i]],0.95)
  sp_info[i,50]=quantile(t_rvc_kist[[i]],0.025)
  sp_info[i,51]=quantile(t_rvc_kist[[i]],0.975)
  sp_info[i,52]=median(t_reef_kist[[i]])
  sp_info[i,53]=quantile(t_reef_kist[[i]],0.1)
  sp_info[i,54]=quantile(t_reef_kist[[i]],0.9)
  sp_info[i,55]=quantile(t_reef_kist[[i]],0.05)
  sp_info[i,56]=quantile(t_reef_kist[[i]],0.95)
  sp_info[i,57]=quantile(t_reef_kist[[i]],0.025)
  sp_info[i,58]=quantile(t_reef_kist[[i]],0.975)
  sp_info[i,59:67]=apply(sp[,grepl('a_hab1',colnames(sp))],2,median)
  sp_info[i,68:72]=apply(sp[,grepl('a_hab2',colnames(sp))],2,median)
  sp_info[i,73:78]=apply(sp[,grepl('a_strat1',colnames(sp))],2,median)
  sp_info[i,79:84]=apply(sp[,grepl('a_strat2',colnames(sp))],2,median)
  
  }

##1. correlations between surveys in year-to-year fluctuations####
r=unlist(r_list)

#find mode of each correlation from posterior
sp_info$mode.r=NA
for(i in 1:length(r_list)){
  d=density(r_list[[i]])
  sp_info$mode.r[i]=d$x[which.max(d$y)]
}

par(mfrow=c(1,2),mar=c(4,1,1,1))
r_m=as.matrix(dplyr::bind_cols(r_list))
colnames(r_m)=seq(1:87)

r_m2=(apply(r_m,1,median)) #median r: 0.377
quantile(apply(r_m,1,median),0.05) #l95 r: 0.265
quantile(apply(r_m,1,median),0.95) #u95: 0.49
hist(r,breaks=30,col='transparent',border='black',xlim=c(-1,1),xlab='',main='',ylab='',yaxt='n',freq=FALSE,lwd=0.5)
par(new=T)
hist(r_m2,breaks=30,col=adjustcolor('black',alpha.f=0.6),border='transparent',xlim=c(-1,1),xlab=expression(rho),main='',ylab='',yaxt='n',freq=FALSE,lwd=0.5)
abline(v=0,lty=5)
mtext(side=3,at=-1,line=-1,substitute(paste(bold('A'))))

par(mar=c(4,1,1,1))
plot(seq(0,88,length.out=12)~seq(-1,1,length.out=12),bty='n',ylab='',xlab=expression(rho),type='n',yaxt='n')
mtext(side=3,at=-1,line=-1,substitute(paste(bold('B'))))
abline(v=0,lty=5)
sp_info2=sp_info %>% arrange(r)
for(i in 1:nrow(sp_info2)){
  
 # lines(rep(i,2)~c(sp_info2$r_u95[i],sp_info2$r_l95[i]),lwd=0.5)
  lines(rep(i,2)~c(sp_info2$r_u90[i],sp_info2$r_l90[i]),lwd=1,col=adjustcolor('black',alpha.f=0.6))
  lines(rep(i,2)~c(sp_info2$r_u80[i],sp_info2$r_l80[i]),lwd=2,col=adjustcolor('black',alpha.f=0.6))
  points(i~sp_info2$r[i],pch=21,bg=adjustcolor('black',alpha.f=0.6),col='white',cex=0.8,lwd=0.5)
  points(i~sp_info2$mode.r[i],pch=24,bg='white',col='black',cex=0.6,lwd=0.5)
#  points(i~sp_info$r[i],pch=21,bg='white',col='black',cex=0.8)
  
}
text(y=83,x=-1,label='CI')
text(y=83,x=-0.8,label='80%')
text(y=83,x=-0.55,label='90%')
lines(rep(80,2)~c(-0.9,-0.75),lwd=2)
lines(rep(80,2)~c(-0.65,-0.5),lwd=1)



##3. by family####
fam_n<- sp_info %>% group_by(Family) %>% summarize(n=n())
fam_n2=subset(fam_n,n>2)

f_info=list()
f_corr=list()
f_r_dist=list()
f_dat=data.frame(grp=fam_n2$Family,med.r=NA,l95=NA,u95=NA,l90=NA,u90=NA,l80=NA,u80=NA,x=seq(1:11),n=fam_n2$n)
for(i in 1:nrow(fam_n2)){
  f_info[[i]]=subset(sp_info,Family==fam_n2$Family[i])
  f_corr[[i]]=r_m[,as.numeric(rownames(f_info[[i]]))]
  f_dat[i,2]=median(apply(f_corr[[i]],1,mean))
  f_dat[i,3]=quantile(apply(f_corr[[i]],1,mean),0.025)
  f_dat[i,4]=quantile(apply(f_corr[[i]],1,mean),0.975)
  f_dat[i,5]=quantile(apply(f_corr[[i]],1,mean),0.05)
  f_dat[i,6]=quantile(apply(f_corr[[i]],1,mean),0.95)
  f_dat[i,7]=quantile(apply(f_corr[[i]],1,mean),0.1)
  f_dat[i,8]=quantile(apply(f_corr[[i]],1,mean),0.9)
  
  f_r_dist[[i]]=density(apply(f_corr[[i]],1,mean),bw=0.02)
  
}

f_dat2= f_dat %>% arrange(med.r)

dev.off()
par(mar=c(4,10,1,2),xpd=T)
plot(seq(0.5,12.5,length.out=12)~seq(-1,1,length.out=12),bty='l',xlab='',ylab='',type='n',yaxt='n',ylim=c(1,12.8))
lines(c(0.5,12.5)~rep(0,2),lty=5)
text(y=12.5,x=-1.2,xpd=T,labels="n",adj=0)
for(i in 1:11){
  f_dist=f_r_dist[[f_dat2$x[i]]]
  f_dist$y=f_dist$y+i
  polygon(f_dist,col=adjustcolor('gray',alpha.f=0.7),border=NA)
  lines(rep(i+0.4,2)~c(f_dat2[i,3],f_dat2[i,4]),lwd=1)
  lines(rep(i+0.4,2)~c(f_dat2[i,5],f_dat2[i,6]),lwd=3)
  lines(rep(i+0.4,2)~c(f_dat2[i,7],f_dat2[i,8]),lwd=5)
#  peak_r=f_dist$x[f_dist$y==max(f_dist$y)]
#  lines(c(min(f_dist$y),max(f_dist$y))~rep(peak_r,2),lwd=1,col='navy')
  points(i+0.4~f_dat2[i,2],pch=21,bg=adjustcolor('black',alpha.f=0.6),col='white',cex=1.5)
  text(y=i+0.8,x=-1.8,xpd=T,labels=f_dat2$grp[i],adj=0)
  text(y=i+0.8,x=-1.2,xpd=T,labels=f_dat2$n[i],adj=0)
}
text(y=12.8,x=-1,label='CI')
text(y=12.8,x=-0.825,label='80%')
text(y=12.8,x=-0.625,label='90%')
text(y=12.8,x=-0.425,label='95%')
lines(rep(12.4,2)~c(-0.9,-0.75),lwd=5)
lines(rep(12.4,2)~c(-0.7,-0.55),lwd=3)
lines(rep(12.4,2)~c(-0.5,-0.35),lwd=1)
mtext(expression(bar(rho)),side=1,line=2.5,at=0)


##4. assess correlation in mean abundance####
cor_kl<- NA
cor_kl_list<- list()
for(q in 1:nrow(m_abund_list[[1]])){
  cor_kl_list[[q]]<- log10(m_abund_list[[1]][q,1:2])
  for(z in 2:length(m_abund_list)){
    cor_kl_list[[q]]<- rbind(cor_kl_list[[q]],log10(m_abund_list[[z]][q,1:2]))
  }   
  cor_kl[q]<- cor.test(cor_kl_list[[q]][,1],cor_kl_list[[q]][,2],method='spearman')$estimate
  }
mean(cor_kl) #mean spearman (for rank-correlation) rho: 0.66
quantile(cor_kl,0.025) #lower quantile
quantile(cor_kl,0.975) #upper quantile

dev.off()

cols=RColorBrewer::brewer.pal(n=11,name='RdYlBu')
br=seq(-1,1,length.out=11)
col_points=cols[findInterval(sp_info$r,br)]

plot(log10(m_abund_reef)~log10(m_abund_rvc),data=sp_info,bty='l',xaxt='n',yaxt='n',xlab='Mean encounter rate (RVC)',ylab='Mean encounter rate (REEF)',type='n',ylim=c(-3,1.5),xlim=c(-3,1.5),cex.lab=1.2)
points(log10(m_abund_reef)~log10(m_abund_rvc),data=sp_info,pch=21,bg=col_points,col='transparent',cex=1.5)
for(i in 1:nrow(sp_info)){
  lines(c(log10(sp_info$m_abund_reef_l90[i]),log10(sp_info$m_abund_reef_u90[i]))~rep(log10(sp_info$m_abund_rvc[i]),2),lwd=1,col=adjustcolor(col_points[i],alpha.f=0.8))
  lines(rep(log10(sp_info$m_abund_reef[i]),2)~c(log10(sp_info$m_abund_rvc_l90[i]),log10(sp_info$m_abund_rvc_u90[i])),lwd=1,col=adjustcolor(col_points[i],alpha.f=0.8))
  
  }

axis(2, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)

axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
median(cor_kl)
quantile(cor_kl,0.025)
quantile(cor_kl,0.975)
text(x=log10(0.005),y=log10(0.001),expression(paste(italic(r),' = 0.66 [0.65 - 0.70]')))

plotdim <- par("plt")
xleft    = plotdim[2] - (plotdim[2] - plotdim[1]) * 0.05
xright   = plotdim[2]  #
ybottom  = plotdim[4] - (plotdim[4] - plotdim[3]) * 0.95  #
ytop     = plotdim[4] - (plotdim[4] - plotdim[3]) * 0.65#

par(
  fig = c(xleft, xright, ybottom, ytop)
  , mar=c(0,0,0,0)
  , new=TRUE
)
x=1
y=br
z=matrix(1:11,nrow=1)
image(x,y,z,col=adjustcolor(cols,alpha.f = 0.9),axes=F,xlab='',ylab='')
axis(2,cex.axis=0.6)
mtext(side=3,expression(tilde(rho)),cex=0.8)

#text(-3,par("usr")[2] + 0.05,'A',cex=1,pos=3,xpd=T,font=2)


#trait-based differences#### Table S2
##3. comparing across species traits####
###color####
col=subset(sp_info,color=='colorful');drab=subset(sp_info,color=='drab')

r_col=r_m[,(colnames(r_m) %in% rownames(col))]
r_drab=r_m[,(colnames(r_m) %in% rownames(drab))]

ncol(r_col)
ncol(r_drab)

median(apply(r_col,1,mean)) #median 
quantile(apply(r_col,1,mean),0.025) #l95 
quantile(apply(r_col,1,mean),0.975) #u95

median(apply(r_drab,1,mean)) #median r
quantile(apply(r_drab,1,mean),0.025) #l95 
quantile(apply(r_drab,1,mean),0.975) #u95

#difference
median(apply(r_col,1,mean)-apply(r_drab,1,mean)) #median diff (state 1:2): -0.125
quantile(apply(r_col,1,mean)-apply(r_drab,1,mean),0.025) #l95: -0.401
quantile(apply(r_col,1,mean)-apply(r_drab,1,mean),0.975) #u95: -0.169


###crypsis####
cryp=subset(sp_info,behavior=='cryptic');consp=subset(sp_info,behavior=='conspicuous')

r_cryp=r_m[,(rownames(cryp))]
r_consp=r_m[,(rownames(consp))]

ncol(r_cryp)
ncol(r_consp)

median(apply(r_consp,1,mean)) #median 
quantile(apply(r_consp,1,mean),0.025) #l95 
quantile(apply(r_consp,1,mean),0.975) #u95

median(apply(r_cryp,1,mean)) #median r: 
quantile(apply(r_cryp,1,mean),0.025) #l95 
quantile(apply(r_cryp,1,mean),0.975) #u95:

#difference
median(apply(r_consp,1,mean)-apply(r_cryp,1,mean)) #median diff (state 1:2): 0.095
quantile(apply(r_consp,1,mean)-apply(r_cryp,1,mean),0.025) #l95: -0.292
quantile(apply(r_consp,1,mean)-apply(r_cryp,1,mean),0.975) #u95: -0.517

###schooling####
sol=subset(sp_info,Schooling==0);schol=subset(sp_info,Schooling==1)

r_sol=r_m[,(rownames(sol))]
r_schol=r_m[,(rownames(schol))]

median(apply(r_sol,1,mean)) #median 
quantile(apply(r_sol,1,mean),0.025) #l95 
quantile(apply(r_sol,1,mean),0.975) #u95

median(apply(r_schol,1,mean)) #median r: 0.195
quantile(apply(r_schol,1,mean),0.025) #l95 r: -0.155
quantile(apply(r_schol,1,mean),0.975) #u95: 0.505

#difference
median(apply(r_sol,1,mean)-apply(r_schol,1,mean)) #median diff (state 1:2): 0.095
quantile(apply(r_sol,1,mean)-apply(r_schol,1,mean),0.025) #l95: -0.292
quantile(apply(r_sol,1,mean)-apply(r_schol,1,mean),0.975) #u95: -0.517

###solitary####
sol01=subset(sp_info,Solitary==0);sol02=subset(sp_info,Solitary==1)

r_sol01=r[,(rownames(sol01))]
r_sol02=r[,(rownames(sol02))]

median(apply(r_sol01,1,median)-apply(r_sol02,1,median)) #median r: 0.377
quantile(apply(r_sol,1,median)-apply(r_schol,1,median),0.025) #l95 r: 0.233
quantile(apply(r_sol,1,median)-apply(r_schol,1,median),0.975) #u95: 0.503

hist(apply(r_sol,1,median)-apply(r_schol,1,median),breaks=30)

median(apply(r_consp,1,median)) #median r: 0.377
quantile(apply(r_drab,1,median),0.025) #l95 r: 0.233
quantile(apply(r_drab,1,median),0.975) #u95: 0.503

###Size####
cor_size=matrix(nrow=nrow(r_m),ncol=3)
for(i in 1:nrow(r_m)){
  cor_size[i,1]=cor.test(r_m[i,],log(sp_info$size))$estimate
  cor_size[i,2]=cor.test(r_m[i,],log(sp_info$size))$conf.int[1]
  cor_size[i,3]=cor.test(r_m[i,],log(sp_info$size))$conf.int[2]
}

median(cor_size[,1])
median(cor_size[,2])
median(cor_size[,3])


###mean abundance####
cor_abund=matrix(nrow=nrow(r_m),ncol=6)
for(i in 1:nrow(r_m)){
  cor_abund[i,1]=cor.test(r_m[i,],cor_kl_list[[i]][,1])$estimate
  cor_abund[i,2]=cor.test(r_m[i,],cor_kl_list[[i]][,1])$conf.int[1]
  cor_abund[i,3]=cor.test(r_m[i,],cor_kl_list[[i]][,1])$conf.int[2]
  cor_abund[i,4]=cor.test(r_m[i,],cor_kl_list[[i]][,2])$estimate
  cor_abund[i,5]=cor.test(r_m[i,],cor_kl_list[[i]][,2])$conf.int[1]
  cor_abund[i,6]=cor.test(r_m[i,],cor_kl_list[[i]][,2])$conf.int[2]
}

median(cor_abund[,1])
median(cor_abund[,2])
median(cor_abund[,3])

median(cor_abund[,4])
median(cor_abund[,5])
median(cor_abund[,6])


#total survey effort
load(here('data','data.Rdata'))

#data overview - supplemental
reef_surv_eff= reef_occs_3403[[1]] %>% group_by(year) %>% summarize(n=n(),n.sites=length(unique(geogr)))
rvc_surv_eff= rvc_occs_3403[[1]] %>% group_by(YEAR) %>% summarize(n=n(),n.sites=length(unique(PRIMARY_SAMPLE_UNIT)))

print(rvc_surv_eff,n=30)
print(reef_surv_eff,n=30)

reef_surv_hab= reef_occs_3403[[1]] %>% group_by(hab_class) %>% summarize(n=n(),n.sites=length(unique(geogr)))
reef_surv_strat= reef_occs_3403[[1]] %>% group_by(stratum) %>% summarize(n=n(),n.sites=length(unique(geogr)))

rvc_surv_hab= rvc_occs_3403[[1]] %>% group_by(HABITAT_CD) %>% summarize(n=n(),n.sites=length(unique(psu_id)))
rvc_surv_strat= rvc_occs_3403[[1]] %>% group_by(STRAT) %>% summarize(n=n(),n.sites=length(unique(psu_id)))


###Variance decomposition####
#For fixed effects (bottom time, depth, etc.) - need variance of effects from data for conditional variance


files_kl_par_2<- list.files(here('outputs','parameter estimates'))

var_frame_reef_kl<- data.frame(sp_info$commonname,prop.site=NA,prop.hab=NA,prop.strat=NA,prop.diver=NA,prop.d.cluster=NA,prop.m.cluster=NA,prop.mth=NA,prop.process=NA,prop.obserror=NA,prop.btime=NA,prop.depth=NA,prop.vis=NA,prop.current=NA,prop.expert=NA)
var_frame_rvc_kl<- data.frame(sp_info$commonname,prop.psu=NA,prop.hab=NA,prop.strat=NA,prop.mth=NA,prop.process=NA,prop.obserror=NA,prop.depth=NA)
for(i in 1:nrow(sp_info)){
  sp<- read.csv(here('outputs','parameter estimates',files_kl_par_2[i]))
  reef_occs<- reef_occs_3403[[i]]
  rvc_occs<- rvc_occs_3403[[i]]
  betas1<-sp[,(gsub('\\..*','',colnames(sp))=='beta1')]
  betas2<-sp[,(gsub('\\..*','',colnames(sp))=='beta2')]
  X1<- matrix(data=c(scale(as.numeric(rvc_occs$DEPTH)),scale(as.numeric(rvc_occs$DEPTH^2))),ncol=2,nrow=nrow(rvc_occs))
  X2<- matrix(data=c(scale(as.numeric(reef_occs$btime)),scale(as.numeric(reef_occs$averagedepth)),scale(as.numeric(reef_occs$averagedepth^2)),scale(as.numeric(reef_occs$visibility)),scale(as.numeric(reef_occs$current)),reef_occs$exp_binary),ncol=6,nrow=nrow(reef_occs))
  var_beta1<- matrix(nrow=nrow(betas1),ncol=ncol(betas1))
  var_beta2<- matrix(nrow=nrow(betas2),ncol=ncol(betas2))
  for(q in 1:ncol(betas1)){
    for(z in 1:nrow(betas1)){ 
      var_beta1[z,q]<- var(betas1[z,q]*X1[,q]) #variance attributed to fixed effects - Nakagawa & Schielzeth 20
    }
  }
  for(q in 1:ncol(betas2)){
      for(z in 1:nrow(betas2)){ 
        var_beta2[z,q]<- var(betas2[z,q]*X2[,q]) #variance attributed to fixed effects - Nakagawa & Schielzeth 20
      }
  }
  var_full1<- data.frame(sp$sd_psu^2,sp$sd_hab1^2,sp$sd_strat1^2,sp$sd_mth1^2,sp$sd_q.1.^2,sp$sd_r.1.^2,var_beta1)
  var_full2<- data.frame(sp$sd_site^2,sp$sd_hab2^2,sp$sd_strat2^2,sp$sd_dv^2,sp$sd_dmy^2,sp$sd_my^2,sp$sd_mth1^2,sp$sd_q.2.^2,sp$sd_r.2.^2,var_beta2)

  var_frame_rvc_kl[i,2]=median(var_full1[,1])  
  var_frame_rvc_kl[i,3]=median(var_full1[,2])  
  var_frame_rvc_kl[i,4]=median(var_full1[,3])  
  var_frame_rvc_kl[i,5]=median(var_full1[,4])
  var_frame_rvc_kl[i,6]=median(var_full1[,5])  
  var_frame_rvc_kl[i,7]=median(var_full1[,6])
  var_frame_rvc_kl[i,8]=median(var_full1[,7]+var_full1[,7]) #depth + depth2
  
  
  var_frame_reef_kl[i,2]=median(var_full2[,1])  
  var_frame_reef_kl[i,3]=median(var_full2[,2])  
  var_frame_reef_kl[i,4]=median(var_full2[,3])  
  var_frame_reef_kl[i,5]=median(var_full2[,4])  
  var_frame_reef_kl[i,6]=median(var_full2[,5])  
  var_frame_reef_kl[i,7]=median(var_full2[,6]) 
  var_frame_reef_kl[i,8]=median(var_full2[,7]) 
  var_frame_reef_kl[i,9]=median(var_full2[,8]) 
  var_frame_reef_kl[i,10]=median(var_full2[,9])
  var_frame_reef_kl[i,11]=median(var_full2[,10])
  var_frame_reef_kl[i,12]=median(var_full2[,11])
  var_frame_reef_kl[i,13]=median(var_full2[,12]+var_full2[,13])#depth + depth2
  var_frame_reef_kl[i,14]=median(var_full2[,14])
  var_frame_reef_kl[i,15]=median(var_full2[,15])

}
write.csv(var_frame_rvc_kl,here('outputs','variance decomp','species_variances_rvc_kl.csv'))
write.csv(var_frame_reef_kl,here('outputs','variance decomp','species_variances_reef_kl.csv'))

var_frame_reef_kl$total.var=rowSums(var_frame_reef_kl[,2:15])
var_frame_reef_kl$spatial=rowSums(var_frame_reef_kl[,2:4])
var_frame_reef_kl$prop.spatial=var_frame_reef_kl$spatial/var_frame_reef_kl$total.var
summary(var_frame_reef_kl$prop.spatial)

var_frame_rvc_kl$total.var=rowSums(var_frame_rvc_kl[,2:8])
var_frame_rvc_kl$spatial=rowSums(var_frame_rvc_kl[,2:4])
var_frame_rvc_kl$prop.spatial=var_frame_rvc_kl$spatial/var_frame_rvc_kl$total.var
summary(var_frame_rvc_kl$prop.spatial)
####