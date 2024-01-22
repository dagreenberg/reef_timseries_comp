#Output analysis
library(here); library(dplyr);library(rlist)

###Functions###
source(here('code','functions.R'))

####data load####
sp_info=read.csv(here('outputs','sp_info.csv'))
dat2=read.csv(here('data','REEF','Caribbean_fish_trait_matrix.csv'))
sp_info$target=dat2$target[match(sp_info$commonname,dat2$commonname)]
sp_info$Family=dat2$Family[match(sp_info$commonname,dat2$commonname)]
pars<- list.files(here('outputs','parameter estimates'))

##check model diagnostics ####
#3 scores: rhat >= 1.05, ess_tail >400, ess_bulk > 400
mod_diag=read.csv(here('outputs','model diagnostics','01_Blue Angelfish_diag.csv'))
for(q in 2:nrow(sp_info)){
  m=read.csv(here('outputs','model diagnostics',paste(sprintf("%02d",q),'_',sp_info$commonname[q],'_diag.csv',sep='')))
  mod_diag=rbind(mod_diag,m)
}
any(mod_diag[,2:4]==T)


m_abund_list<- list() #mean abundance
r_list<- list() #temp. correlation
g_rvc_list<- list();t_rvc_list=list()
g_reef_list<- list();t_reef_list=list()

#Overview####
sp_info$m_abund_rvc=NA;sp_info$m_abund_rvc_l90=NA;sp_info$m_abund_rvc_u90=NA;sp_info$m_abund_rvc_l95=NA;sp_info$m_abund_rvc_u95=NA
sp_info$m_abund_reef=NA;sp_info$m_abund_reef_l90=NA;sp_info$m_abund_reef_u90=NA;sp_info$m_abund_reef_l95=NA;sp_info$m_abund_reef_u95=NA
sp_info$r=NA;sp_info$r.prob=NA;sp_info$r_l80=NA;sp_info$r_u80=NA;sp_info$r_l90=NA;sp_info$r_u90=NA;sp_info$r_l95=NA;sp_info$r_u95=NA
sp_info$trend_rvc=NA;sp_info$trend_rvc_l80=NA;sp_info$trend_rvc_u80=NA;sp_info$trend_rvc_l90=NA;sp_info$trend_rvc_u90=NA;sp_info$trend_rvc_l95=NA;sp_info$trend_rvc_u95=NA
sp_info$trend_reef=NA;sp_info$trend_reef_l80=NA;sp_info$trend_reef_u80=NA;sp_info$trend_reef_l90=NA;sp_info$trend_reef_u90=NA;sp_info$trend_reef_l95=NA;sp_info$trend_reef_u95=NA
sp_info$a_hab_1.1=NA;sp_info$a_hab_1.2=NA;sp_info$a_hab_1.3=NA;sp_info$a_hab_1.4=NA;sp_info$a_hab_1.5=NA;sp_info$a_hab_1.6=NA;sp_info$a_hab_1.7=NA;sp_info$a_hab_1.8=NA;sp_info$a_hab_1.9=NA
sp_info$a_hab_2.1=NA;sp_info$a_hab_2.2=NA;sp_info$a_hab_2.3=NA;sp_info$a_hab_2.4=NA;sp_info$a_hab_2.5=NA;
sp_info$a_strat_1.1=NA;sp_info$a_strat_1.2=NA;sp_info$a_strat_1.3=NA;sp_info$a_strat_1.4=NA;sp_info$a_strat_1.5=NA;sp_info$a_strat_1.6=NA;sp_info$a_strat_1.7=NA
sp_info$a_strat_2.1=NA;sp_info$a_strat_2.2=NA;sp_info$a_strat_2.3=NA;sp_info$a_strat_2.4=NA;sp_info$a_strat_2.5=NA;sp_info$a_strat_2.6=NA
sp_info$diff.abund=NA;sp_info$diff.abund_l90=NA;sp_info$diff.abund_u90=NA;sp_info$diff.abund_l95=NA;sp_info$diff.abund_u95=NA
sp_info$phi=NA;sp_info$sd_r.rvc=NA
for(i in 1:nrow(sp_info)){
  sp<- read.csv(here('outputs','parameter estimates',pars[i]))
  x_mat=sp[,grepl('x',colnames(sp))]
  x1=x_mat[,1:26]
  x2=x_mat[,27:ncol(x_mat)]
  c=sp[,grepl('cut',colnames(sp))]
 #mean abundance among all years:
  m_abund_rvc<- exp(apply(x1,1,mean))
  m_abund_reef<- exp(apply(log(mean_ord_to_n(x=x2,c=c)),1,mean))
  m_abund_list[[i]]<- data.frame(mean.rvc.kl=sort(m_abund_rvc),mean.reef.kl=sort(m_abund_reef),sp=sp_info$commonname[i],mean.diff=sort(m_abund_reef-m_abund_rvc))
  #population trend: 
  g_rvc_list[[i]]<- matrix(data=NA,nrow=4800,ncol=ncol(x1)-1)
  t_rvc_list[[i]]<- rep(0,4800)
  for(t in 1:nrow(x1)){
    for(z in 1:(ncol(x1)-1)){
      g_rvc_list[[i]][t,z]=x1[t,z+1]-x1[t,z]
    }
    t_rvc_list[[i]][t]<- (log(gm_mean(exp(g_rvc_list[[i]][t,]))))
  }
  
  x_abund_reef<- mean_ord_to_n(x=x2,c=c)
  g_reef_list[[i]]<- matrix(data=NA,nrow=4800,ncol=ncol(x2)-1)
  t_reef_list[[i]]<- rep(0,4800)
  for(t in 1:nrow(x2)){
    for(z in 1:(ncol(x2)-1)){
      g_reef_list[[i]][t,z]=log(x_abund_reef[t,z+1])-log(x_abund_reef[t,z])
    }
    t_reef_list[[i]][t]<- log(gm_mean(exp(g_reef_list[[i]][t,])))
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
  sp_info[i,45]=median(t_rvc_list[[i]])
  sp_info[i,46]=quantile(t_rvc_list[[i]],0.1)
  sp_info[i,47]=quantile(t_rvc_list[[i]],0.9)
  sp_info[i,48]=quantile(t_rvc_list[[i]],0.05)
  sp_info[i,49]=quantile(t_rvc_list[[i]],0.95)
  sp_info[i,50]=quantile(t_rvc_list[[i]],0.025)
  sp_info[i,51]=quantile(t_rvc_list[[i]],0.975)
  sp_info[i,52]=median(t_reef_list[[i]])
  sp_info[i,53]=quantile(t_reef_list[[i]],0.1)
  sp_info[i,54]=quantile(t_reef_list[[i]],0.9)
  sp_info[i,55]=quantile(t_reef_list[[i]],0.05)
  sp_info[i,56]=quantile(t_reef_list[[i]],0.95)
  sp_info[i,57]=quantile(t_reef_list[[i]],0.025)
  sp_info[i,58]=quantile(t_reef_list[[i]],0.975)
  sp_info[i,59:67]=apply(sp[,grepl('a_hab1',colnames(sp))],2,median)
  sp_info[i,68:72]=apply(sp[,grepl('a_hab2',colnames(sp))],2,median)
  sp_info[i,73:79]=apply(sp[,grepl('a_strat1',colnames(sp))],2,median)
  sp_info[i,80:85]=apply(sp[,grepl('a_strat2',colnames(sp))],2,median)
  sp_info[i,86]=median(log10(m_abund_reef)-log10(m_abund_rvc))
  sp_info[i,87]=quantile(log10(m_abund_reef)-log10(m_abund_rvc),0.05)
  sp_info[i,88]=quantile(log10(m_abund_reef)-log10(m_abund_rvc),0.95)
  sp_info[i,89]=quantile(log10(m_abund_reef)-log10(m_abund_rvc),0.025)
  sp_info[i,90]=quantile(log10(m_abund_reef)-log10(m_abund_rvc),0.975)
  sp_info[i,91]=median(sp[,grepl('phi',colnames(sp))])
  sp_info[i,92]=median(sp[,grepl('sd_r.1.',colnames(sp))])
  }



#overall correlation among species####
r=unlist(r_list)

#find mode of each correlation from posterior
sp_info$mode.r=NA
for(i in 1:length(r_list)){
  d=density(r_list[[i]])
  sp_info$mode.r[i]=d$x[which.max(d$y)]
}
#summary of modes for rho
summary(sp_info$mode.r)
quantile(sp_info$mode.r,seq(0,1,by=0.1))

#summary of medians for rho
summary(sp_info$r)
quantile(sp_info$r,seq(0,1,by=0.1))

r_m=as.matrix(dplyr::bind_cols(r_list))
colnames(r_m)=seq(1:87)

## taxonomic differences ####
fam_cols2=c('#a50026',
  '#d73027',
  '#f46d43',
  '#fdae61',
  '#fee090',
  '#ffffbf',
  '#f6e8c3',
  '#e0f3f8',
  '#abd9e9',
  '#74add1',
  '#4575b4',
  '#313695')

#keep all families with 3 sp or more
fam_n<- sp_info %>% group_by(Family) %>% summarize(n=n())
fam_n2=subset(fam_n,n>2)

f_info=list()
f_corr=list()
f_r_dist=list()
f_dat=data.frame(grp=fam_n2$Family,med.r=NA,l95=NA,u95=NA,l90=NA,u90=NA,l80=NA,u80=NA,x=seq(1:12),n=fam_n2$n)
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

f_dat$cols=rev(fam_cols2)
f_dat2= f_dat %>% arrange(med.r)
f_dat2$cols=fam_cols2

dev.off()

sp_info$f.cols=f_dat2$cols[match(sp_info$Family,f_dat2$grp)]
sp_info$f.cols=ifelse(is.na(sp_info$f.cols),'darkgray',sp_info$f.cols)

par(mfrow=c(1,2),mar=c(4,1,1,0))
plot(seq(0,88,length.out=12)~seq(-1,1,length.out=12),bty='n',ylab='',xlab='',type='n',yaxt='n')
mtext(side=3,at=-1,line=-0.3,substitute(paste(bold('A'))))
mtext(expression(rho),side=1,line=2.5,at=0)
lines(c(0,88)~rep(0,2),lty=5)
sp_info2=sp_info %>% arrange(r)
for(i in 1:nrow(sp_info2)){
 # lines(rep(i,2)~c(sp_info2$r_u95[i],sp_info2$r_l95[i]),lwd=0.5)
  lines(rep(i,2)~c(sp_info2$r_u90[i],sp_info2$r_l90[i]),lwd=1,col=adjustcolor(sp_info2$f.cols[i],alpha.f=0.6))
  lines(rep(i,2)~c(sp_info2$r_u80[i],sp_info2$r_l80[i]),lwd=2,col=adjustcolor(sp_info2$f.cols[i],alpha.f=0.6))
  points(i~sp_info2$r[i],pch=21,bg=sp_info2$f.cols[i],col='black',cex=0.8,lwd=0.5)
  points(i~sp_info2$mode.r[i],pch=24,bg=sp_info2$f.cols[i],col='black',cex=0.7,lwd=0.5)
  #points(i~sp_info$r[i],pch=21,bg='white',col='black',cex=0.8)
  
}
text(y=par('usr')[4]*0.95,x=-1,label='CI',cex=0.8)
text(y=par('usr')[4]*0.95,x=-0.8,label='80%',cex=0.8)
text(y=par('usr')[4]*0.95,x=-0.6,label='90%',cex=0.8)
lines(rep(par('usr')[4]*0.93,2)~c(-0.9,-0.75),lwd=2)
lines(rep(par('usr')[4]*0.93,2)~c(-0.7,-0.55),lwd=1)

par(mar=c(4,1,1,1),xpd=T)
plot(seq(0.5,13.5,length.out=12)~seq(-1,1,length.out=12),bty='n',xlab='',ylab='',type='n',yaxt='n',ylim=c(1,13.8))
mtext(side=3,at=-1,line=-0.3,substitute(paste(bold('B'))))
lines(c(0.5,13.5)~rep(0,2),lty=5)
text(y=13.5,x=1.05,xpd=T,labels=substitute(paste(bold('n'))),adj=0,cex=0.8)
for(i in 12:1){
  f_dist=f_r_dist[[f_dat2$x[i]]]
  f_dist$y=f_dist$y+i
  polygon(f_dist,col=adjustcolor(f_dat2$cols[i],alpha.f=0.6),border=NA)
  lines(rep(i+0.4,2)~c(f_dat2[i,3],f_dat2[i,4]),lwd=1)
  lines(rep(i+0.4,2)~c(f_dat2[i,5],f_dat2[i,6]),lwd=3)
#  lines(rep(i+0.4,2)~c(f_dat2[i,7],f_dat2[i,8]),lwd=5)
#  peak_r=f_dist$x[f_dist$y==max(f_dist$y)]
#  lines(c(min(f_dist$y),max(f_dist$y))~rep(peak_r,2),lwd=1,col='navy')

  text(y=i+0.5,x=-1.2,xpd=T,labels=f_dat2$grp[i],adj=0,cex=0.8)
  text(y=i+0.5,x=1.05,xpd=T,labels=f_dat2$n[i],adj=0,cex=0.8)
}
for(i in 12:1){
  points(i+0.4~f_dat2[i,2],pch=21,bg='black',col='white',cex=1.5)
}
text(y=par('usr')[4]*0.95,x=-1,label='CI',cex=0.8)
text(y=par('usr')[4]*0.95,x=-0.8,label='80%',cex=0.8)
text(y=par('usr')[4]*0.95,x=-0.6,label='90%',cex=0.8)
#text(y=par('usr')[4]*0.93,x=-0.4,label='95%',cex=0.8)
lines(rep(par('usr')[4]*0.93,2)~c(-0.9,-0.75),lwd=3)
lines(rep(par('usr')[4]*0.93,2)~c(-0.7,-0.55),lwd=1)
#lines(rep(par('usr')[4]*0.91,2)~c(-0.5,-0.35),lwd=1)
mtext(expression(bar(rho)),side=1,line=2.5,at=0)


##survey mean abundance correlation####
cor_kl<- NA
cor_kl_l95<-NA
cor_kl_u95<-NA
cor_kl_list<- list()
for(q in 1:nrow(m_abund_list[[1]])){
  cor_kl_list[[q]]<- log10(m_abund_list[[1]][q,1:2])
  for(z in 2:length(m_abund_list)){
    cor_kl_list[[q]]<- rbind(cor_kl_list[[q]],log10(m_abund_list[[z]][q,1:2]))
  }   
  cor_kl[q]<- cor.test(cor_kl_list[[q]][,1],cor_kl_list[[q]][,2],method='pearson')$estimate
  cor_kl_l95[q]<- cor.test(cor_kl_list[[q]][,1],cor_kl_list[[q]][,2],method='pearson')$conf.int[1]
  cor_kl_u95[q]<- cor.test(cor_kl_list[[q]][,1],cor_kl_list[[q]][,2],method='pearson')$conf.int[2]
  }
median(cor_kl) #mean spearman (for rank-correlation) rho: 0.67
#fully uncertainty:
quantile(cor_kl_l95,0.025) #lower quantile across distribution
quantile(cor_kl_u95,0.975) #upper quantile across distribution

dev.off()

cols=RColorBrewer::brewer.pal(n=11,name='RdYlBu')
br=seq(-1,1,length.out=11)
col_points=cols[findInterval(sp_info$mode.r,br)]

plot(log10(m_abund_reef)~log10(m_abund_rvc),data=sp_info,bty='l',xaxt='n',yaxt='n',xlab='Mean encounter rate (RVC)',ylab='Mean encounter rate (REEF)',type='n',ylim=c(-3,1.5),xlim=c(-3,1.5),cex.lab=1.2)
abline(c(0,1))
for(i in 1:nrow(sp_info)){
  lines(c(log10(sp_info$m_abund_reef_l90[i]),log10(sp_info$m_abund_reef_u90[i]))~rep(log10(sp_info$m_abund_rvc[i]),2),lwd=1,col=adjustcolor(col_points[i],alpha.f=0.3))
  lines(rep(log10(sp_info$m_abund_reef[i]),2)~c(log10(sp_info$m_abund_rvc_l90[i]),log10(sp_info$m_abund_rvc_u90[i])),lwd=1,col=adjustcolor(col_points[i],alpha.f=0.3))
  
}
points(log10(m_abund_reef)~log10(m_abund_rvc),data=sp_info,pch=21,bg=col_points,col='transparent',cex=1.25)

#points(seq(log10(0.3),log10(3),length.out=5)~log10(rep(35,5)),pch=21,cex=quantile(log(sp_info$size)*0.5))
#text(x=log10(35),y=log10(6),'Total length (cm)',cex=0.8,xpd=T)
#text(x=rep(log10(20),5),y=seq(log10(0.3),log10(3),length.out=5),exp(quantile(log(sp_info$size))),cex=0.6)


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
text(x=log10(0.008),y=log10(0.001),expression(paste(italic(r),' = 0.67 [95% CI: 0.57 - 0.79]')))

plotdim <- par("plt")
xleft    = plotdim[2] - (plotdim[2] - plotdim[1]) * 0.05
xright   = plotdim[2]  #
ybottom  = plotdim[4] - (plotdim[4] - plotdim[3]) *1.05 #
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
mtext(side=3,expression(rho),cex=0.8)

#text(-3,par("usr")[2] + 0.05,'A',cex=1,pos=3,xpd=T,font=2)


##traits####

#### color - 2 categories####
d.drab=data.frame(level=c('colourful','drab','col-drab'),med=NA,l90=NA,u90=NA,l80=NA,u80=NA)

drab=subset(sp_info,color=='drab');colf=subset(sp_info,color=='colorful')
r_drab=r_m[,(rownames(drab))]
r_colf=r_m[,(rownames(colf))]

#medians
d.drab[1,2]=median(apply(r_colf,1,mean)) #median 
d.drab[2,2]=median(apply(r_drab,1,mean)) #median 
d.drab[3,2]=median(apply(r_colf,1,mean)-apply(r_drab,1,mean)) #median 

#l90
d.drab[1,3]=quantile(apply(r_colf,1,mean),0.05) #l90 
d.drab[2,3]=quantile(apply(r_drab,1,mean),0.05) #l90
d.drab[3,3]=quantile(apply(r_colf,1,mean)-apply(r_drab,1,mean),0.05) #l90 

#u90
d.drab[1,4]=quantile(apply(r_colf,1,mean),0.95) #u90 
d.drab[2,4]=quantile(apply(r_drab,1,mean),0.95) #u90
d.drab[3,4]=quantile(apply(r_colf,1,mean)-apply(r_drab,1,mean),0.95) #u90 

#l80
d.drab[1,5]=quantile(apply(r_colf,1,mean),0.1) #l80 
d.drab[2,5]=quantile(apply(r_drab,1,mean),0.1) #l80
d.drab[3,5]=quantile(apply(r_colf,1,mean)-apply(r_drab,1,mean),0.1) #l80 

#u80
d.drab[1,6]=quantile(apply(r_colf,1,mean),0.9) #u80 
d.drab[2,6]=quantile(apply(r_drab,1,mean),0.9) #u80
d.drab[3,6]=quantile(apply(r_colf,1,mean)-apply(r_drab,1,mean),0.9) #u80 

print(d.drab)


#### crypsis - 2 categories####
d.cryp=data.frame(level=c('conspicuous','cryptic','consp-cryptic'),med=NA,l90=NA,u90=NA,l80=NA,u80=NA)

cryp=subset(sp_info,behavior=='cryptic');consp=subset(sp_info,behavior=='conspicuous')
r_cryp=r_m[,(rownames(cryp))]
r_consp=r_m[,(rownames(consp))]

#medians
d.cryp[1,2]=median(apply(r_consp,1,mean)) #median 
d.cryp[2,2]=median(apply(r_cryp,1,mean)) #median 
d.cryp[3,2]=median(apply(r_consp,1,mean)-apply(r_cryp,1,mean)) #median 

#l90
d.cryp[1,3]=quantile(apply(r_consp,1,mean),0.05) #l90 
d.cryp[2,3]=quantile(apply(r_cryp,1,mean),0.05) #l90
d.cryp[3,3]=quantile(apply(r_consp,1,mean)-apply(r_cryp,1,mean),0.05) #l90 

#u90
d.cryp[1,4]=quantile(apply(r_consp,1,mean),0.95) #u90 
d.cryp[2,4]=quantile(apply(r_cryp,1,mean),0.95) #u90
d.cryp[3,4]=quantile(apply(r_consp,1,mean)-apply(r_cryp,1,mean),0.95) #u90 

#l80
d.cryp[1,5]=quantile(apply(r_consp,1,mean),0.1) #l80 
d.cryp[2,5]=quantile(apply(r_cryp,1,mean),0.1) #l80
d.cryp[3,5]=quantile(apply(r_consp,1,mean)-apply(r_cryp,1,mean),0.1) #l80 

#u80
d.cryp[1,6]=quantile(apply(r_consp,1,mean),0.9) #u80 
d.cryp[2,6]=quantile(apply(r_cryp,1,mean),0.9) #u80
d.cryp[3,6]=quantile(apply(r_consp,1,mean)-apply(r_cryp,1,mean),0.9) #u80 

#credible intervals for means & mean differences:
print(d.cryp)


#### aggregation - 3 categories####
d.schol=data.frame(level=c('solitary','shoal','school','solitary-shoal','solitary-school','shoal-school'),med=NA,l90=NA,u90=NA,l80=NA,u80=NA)
sol=subset(sp_info,Schooling==0&Solitary==1);shol=subset(sp_info,Schooling==0&Solitary==0);schol=subset(sp_info,Schooling==1)
r_sol=r_m[,(rownames(sol))]
r_shol=r_m[,(rownames(shol))]
r_schol=r_m[,(rownames(schol))]

#median of group means/mean diff.
d.schol[1,2]=median(apply(r_sol,1,mean))
d.schol[2,2]=median(apply(r_shol,1,mean))
d.schol[3,2]=median(apply(r_schol,1,mean))
d.schol[4,2]=median(apply(r_sol,1,mean)-apply(r_shol,1,mean))
d.schol[5,2]=median(apply(r_sol,1,mean)-apply(r_schol,1,mean))
d.schol[6,2]=median(apply(r_shol,1,mean)-apply(r_schol,1,mean))
#90% cis:
d.schol[1,3]=quantile(apply(r_sol,1,mean),0.05)
d.schol[2,3]=quantile(apply(r_shol,1,mean),0.05)
d.schol[3,3]=quantile(apply(r_schol,1,mean),0.05)
d.schol[4,3]=quantile(apply(r_sol,1,mean)-apply(r_shol,1,mean),0.05)
d.schol[5,3]=quantile(apply(r_sol,1,mean)-apply(r_schol,1,mean),0.05)
d.schol[6,3]=quantile(apply(r_shol,1,mean)-apply(r_schol,1,mean),0.05)

d.schol[1,4]=quantile(apply(r_sol,1,mean),0.95)
d.schol[2,4]=quantile(apply(r_shol,1,mean),0.95)
d.schol[3,4]=quantile(apply(r_schol,1,mean),0.95)
d.schol[4,4]=quantile(apply(r_sol,1,mean)-apply(r_shol,1,mean),0.95)
d.schol[5,4]=quantile(apply(r_sol,1,mean)-apply(r_schol,1,mean),0.95)
d.schol[6,4]=quantile(apply(r_shol,1,mean)-apply(r_schol,1,mean),0.95)

#80% cis:
d.schol[1,5]=quantile(apply(r_sol,1,mean),0.1)
d.schol[2,5]=quantile(apply(r_shol,1,mean),0.1)
d.schol[3,5]=quantile(apply(r_schol,1,mean),0.1)
d.schol[4,5]=quantile(apply(r_sol,1,mean)-apply(r_shol,1,mean),0.1)
d.schol[5,5]=quantile(apply(r_sol,1,mean)-apply(r_schol,1,mean),0.1)
d.schol[6,5]=quantile(apply(r_shol,1,mean)-apply(r_schol,1,mean),0.1)

d.schol[1,6]=quantile(apply(r_sol,1,mean),0.9)
d.schol[2,6]=quantile(apply(r_shol,1,mean),0.9)
d.schol[3,6]=quantile(apply(r_schol,1,mean),0.9)
d.schol[4,6]=quantile(apply(r_sol,1,mean)-apply(r_shol,1,mean),0.9)
d.schol[5,6]=quantile(apply(r_sol,1,mean)-apply(r_schol,1,mean),0.9)
d.schol[6,6]=quantile(apply(r_shol,1,mean)-apply(r_schol,1,mean),0.9)

print(d.schol)

###size correlation####
cor_size=matrix(nrow=nrow(r_m),ncol=3)
for(i in 1:nrow(r_m)){
  cor_size[i,1]=cor.test(r_m[i,],log10(sp_info$size))$estimate
  cor_size[i,2]=cor.test(r_m[i,],log10(sp_info$size))$conf.int[1]
  cor_size[i,3]=cor.test(r_m[i,],log10(sp_info$size))$conf.int[2]

}

median(cor_size[,1])
quantile(cor_size[,2],0.025)
quantile(cor_size[,3],0.975)


quantile(sp_info$size,seq(0,1,by=0.2))

d.size=data.frame(level=c('sm','med','lg','sm-med','sm-lg','med-lg'),med=NA,l90=NA,u90=NA,l80=NA,u80=NA)
sm=subset(sp_info,size<=16.4);med=subset(sp_info,size>16.4&size<72.4);lg=subset(sp_info,size>=72.4)
r_sm=r_m[,(rownames(sm))]
r_med=r_m[,(rownames(med))]
r_lg=r_m[,(rownames(lg))]

#median of group means/mean diff.
d.size[1,2]=median(apply(r_sm,1,mean))
d.size[2,2]=median(apply(r_med,1,mean))
d.size[3,2]=median(apply(r_lg,1,mean))
d.size[4,2]=median(apply(r_sm,1,mean)-apply(r_med,1,mean))
d.size[5,2]=median(apply(r_sm,1,mean)-apply(r_lg,1,mean))
d.size[6,2]=median(apply(r_med,1,mean)-apply(r_lg,1,mean))
#90% cis:
d.size[1,3]=quantile(apply(r_sm,1,mean),0.05)
d.size[2,3]=quantile(apply(r_med,1,mean),0.05)
d.size[3,3]=quantile(apply(r_lg,1,mean),0.05)
d.size[4,3]=quantile(apply(r_sm,1,mean)-apply(r_med,1,mean),0.05)
d.size[5,3]=quantile(apply(r_sm,1,mean)-apply(r_lg,1,mean),0.05)
d.size[6,3]=quantile(apply(r_med,1,mean)-apply(r_lg,1,mean),0.05)

d.size[1,4]=quantile(apply(r_sm,1,mean),0.95)
d.size[2,4]=quantile(apply(r_med,1,mean),0.95)
d.size[3,4]=quantile(apply(r_lg,1,mean),0.95)
d.size[4,4]=quantile(apply(r_sm,1,mean)-apply(r_med,1,mean),0.95)
d.size[5,4]=quantile(apply(r_sm,1,mean)-apply(r_lg,1,mean),0.95)
d.size[6,4]=quantile(apply(r_med,1,mean)-apply(r_lg,1,mean),0.95)

#80% cis:
d.size[1,5]=quantile(apply(r_sm,1,mean),0.1)
d.size[2,5]=quantile(apply(r_med,1,mean),0.1)
d.size[3,5]=quantile(apply(r_lg,1,mean),0.1)
d.size[4,5]=quantile(apply(r_sm,1,mean)-apply(r_med,1,mean),0.1)
d.size[5,5]=quantile(apply(r_sm,1,mean)-apply(r_lg,1,mean),0.1)
d.size[6,5]=quantile(apply(r_med,1,mean)-apply(r_lg,1,mean),0.1)

d.size[1,6]=quantile(apply(r_sm,1,mean),0.9)
d.size[2,6]=quantile(apply(r_med,1,mean),0.9)
d.size[3,6]=quantile(apply(r_lg,1,mean),0.9)
d.size[4,6]=quantile(apply(r_sm,1,mean)-apply(r_med,1,mean),0.9)
d.size[5,6]=quantile(apply(r_sm,1,mean)-apply(r_lg,1,mean),0.9)
d.size[6,6]=quantile(apply(r_med,1,mean)-apply(r_lg,1,mean),0.9)

print(d.size)


#figure 4: trait plots####
dev.off()
par(mar=c(4,1,0,1))
plot(seq(0,5,length.out=12)~seq(-0.3,0.8,length.out=12),bty='n',xlab='',ylab='',type='n',yaxt='n',ylim=c(0,5))

dens_col=density(apply(r_colf,1,mean),bw=0.02)
dens_drb=density(apply(r_drab,1,mean),bw=0.02)
dens_col$y=dens_col$y/max(dens_col$y)+3.6
dens_drb$y=dens_drb$y/max(dens_drb$y)+3.6
polygon(dens_drb,col=adjustcolor('#993404',alpha.f=0.6),border=NA)
polygon(dens_col,col=adjustcolor('#fe9929',alpha.f=0.6),border=NA)
text(x=dens_col$x[which.max(dens_col$y)]*0.62,y=dens_col$y[which.max(dens_col$y)],'Colorful',col='#fe9929')
text(x=dens_drb$x[which.max(dens_drb$y)]*1.2,y=dens_drb$y[which.max(dens_drb$y)],'Drab',col='#993404')
lines(rep(3.75,2)~c(d.drab[1,3],d.drab[1,4]),lwd=1)
lines(rep(3.85,2)~c(d.drab[2,3],d.drab[2,4]),lwd=1)
lines(rep(3.75,2)~c(d.drab[1,5],d.drab[1,6]),lwd=3)
lines(rep(3.85,2)~c(d.drab[2,5],d.drab[2,6]),lwd=3)
points(y=c(3.75,3.85),x=d.drab[1:2,2],pch=21,bg=c('#fe9929','#993404'),col='white',cex=1.5)

text(y=4,x=-0.2,label='Coloration',cex=0.8)

dens_cryp=density(apply(r_cryp,1,mean),bw=0.02)
dens_consp=density(apply(r_consp,1,mean),bw=0.02)
dens_consp$y=dens_consp$y/max(dens_consp$y)+2.4
dens_cryp$y=dens_cryp$y/max(dens_cryp$y)+2.4
polygon(dens_cryp,col=adjustcolor('#0868ac',alpha.f=0.6),border=NA)
polygon(dens_consp,col=adjustcolor('#7bccc4',alpha.f=0.6),border=NA)
text(x=dens_consp$x[which.max(dens_consp$y)]*1.4,y=dens_consp$y[which.max(dens_consp$y)],'Conspicuous',col='#7bccc4')
text(x=dens_cryp$x[which.max(dens_cryp$y)]*0.5,y=dens_cryp$y[which.max(dens_cryp$y)],'Cryptic',col='#0868ac')
lines(rep(2.55,2)~c(d.cryp[1,3],d.cryp[1,4]),lwd=1)
lines(rep(2.65,2)~c(d.cryp[2,3],d.cryp[2,4]),lwd=1)
lines(rep(2.55,2)~c(d.cryp[1,5],d.cryp[1,6]),lwd=3)
lines(rep(2.65,2)~c(d.cryp[2,5],d.cryp[2,6]),lwd=3)
points(y=c(2.55,2.65),x=d.cryp[1:2,2],pch=21,bg=c('#7bccc4','#0868ac'),col='white',cex=1.5)

text(y=2.8,x=-0.2,label='Gregariousness',cex=0.8)


dens_sol=density(apply(r_sol,1,mean),bw=0.02)
dens_shol=density(apply(r_shol,1,mean),bw=0.02)
dens_schol=density(apply(r_schol,1,mean),bw=0.02)
dens_sol$y=dens_sol$y/max(dens_sol$y)+1.2
dens_shol$y=dens_shol$y/max(dens_shol$y)+1.2
dens_schol$y=dens_schol$y/max(dens_schol$y)+1.2
polygon(dens_sol,col=adjustcolor('#016c59',alpha.f=0.6),border=NA)

polygon(dens_shol,col=adjustcolor('#1c9099',alpha.f=0.6),border=NA)
polygon(dens_schol,col=adjustcolor('#67a9cf',alpha.f=0.6),border=NA)
text(x=dens_schol$x[which.max(dens_schol$y)]*0.25,y=dens_schol$y[which.max(dens_schol$y)]*1.04,'Schooling',col='#67a9cf')
text(x=dens_shol$x[which.max(dens_shol$y)]*1.05,y=dens_shol$y[which.max(dens_shol$y)]*1.04,'Shoaling',col='#1c9099')
text(x=dens_sol$x[which.max(dens_sol$y)]*1.25,y=dens_sol$y[which.max(dens_sol$y)]*1.04,'Solitary',col='#016c59')
lines(rep(1.3,2)~c(d.schol[1,3],d.schol[1,4]),lwd=1)
lines(rep(1.4,2)~c(d.schol[2,3],d.schol[2,4]),lwd=1)
lines(rep(1.5,2)~c(d.schol[3,3],d.schol[3,4]),lwd=1)
lines(rep(1.3,2)~c(d.schol[1,5],d.schol[1,6]),lwd=3)
lines(rep(1.4,2)~c(d.schol[2,5],d.schol[2,6]),lwd=3)
lines(rep(1.5,2)~c(d.schol[3,5],d.schol[3,6]),lwd=3)
points(y=c(1.3,1.4,1.5),x=d.schol[1:3,2],pch=21,bg=c('#016c59','#1c9099','#67a9cf'),col='white',cex=1.5)

text(y=1.6,x=-0.2,label='Aggregation',cex=0.8)

dens_sm=density(apply(r_sm,1,mean),bw=0.02)
dens_med=density(apply(r_med,1,mean),bw=0.02)
dens_lg=density(apply(r_lg,1,mean),bw=0.02)
dens_sm$y=dens_sm$y/max(dens_sm$y)
dens_med$y=dens_med$y/max(dens_med$y)
dens_lg$y=dens_lg$y/max(dens_lg$y)
polygon(dens_lg,col=adjustcolor('#253494',alpha.f=0.6),border=NA)

polygon(dens_med,col=adjustcolor('#2c7fb8',alpha.f=0.6),border=NA)
polygon(dens_sm,col=adjustcolor('#41b6c4',alpha.f=0.6),border=NA)
text(x=dens_sm$x[which.max(dens_sm$y)]*0.71,y=dens_sm$y[which.max(dens_sm$y)]*1.08,'Small',col='#41b6c4')
text(x=dens_med$x[which.max(dens_med$y)]*1.2,y=dens_med$y[which.max(dens_med$y)]*1.08,'Medium',col='#2c7fb8')
text(x=dens_lg$x[which.max(dens_lg$y)]*1.25,y=dens_lg$y[which.max(dens_lg$y)]*1.06,'Large',col='#253494')
lines(rep(0.1,2)~c(d.size[1,3],d.size[1,4]),lwd=1)
lines(rep(0.2,2)~c(d.size[2,3],d.size[2,4]),lwd=1)
lines(rep(0.3,2)~c(d.size[3,3],d.size[3,4]),lwd=1)
lines(rep(0.1,2)~c(d.size[1,5],d.size[1,6]),lwd=3)
lines(rep(0.2,2)~c(d.size[2,5],d.size[2,6]),lwd=3)
lines(rep(0.3,2)~c(d.size[3,5],d.size[3,6]),lwd=3)
points(y=c(0.1,0.2,0.3),x=d.size[1:3,2],pch=21,bg=c('#41b6c4','#2c7fb8','#253494'),col='white',cex=1.25)
mtext(expression(bar(rho)),side=1,line=2.5)
text(y=0.4,x=-0.2,label='Body Size',cex=0.8)

text(y=par('usr')[4]*0.9,x=-0.25,label='CI',cex=0.8)
text(y=par('usr')[4]*0.9,x=-0.18,label='80%',cex=0.8)
text(y=par('usr')[4]*0.9,x=-0.1,label='90%',cex=0.8)
lines(rep(par('usr')[4]*0.87,2)~c(-0.21,-0.15),lwd=3)
lines(rep(par('usr')[4]*0.87,2)~c(-0.13,-0.07),lwd=1)

###mean abundance trends####

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


dev.off()
par(mfrow=c(3,1))
plot(r~log10(m_abund_rvc),data=sp_info,bty='l',xaxt='n',xlab='Mean encounter rate (RVC)',ylab=expression(rho),type='n',ylim=c(-1,1),xlim=c(-3,1.5),cex.lab=1.2)
for(i in 1:nrow(sp_info)){
  lines(c(sp_info$r_l90[i],sp_info$r_u90[i])~rep(log10(sp_info$m_abund_rvc[i]),2),lwd=1,col=adjustcolor(col_points[i],alpha.f=0.3))
  lines(rep(sp_info$r[i],2)~c(log10(sp_info$m_abund_rvc_l90[i]),log10(sp_info$m_abund_rvc_u90[i])),lwd=1,col=adjustcolor(col_points[i],alpha.f=0.3))
}
points(r~log10(m_abund_rvc),data=sp_info,pch=21,bg=col_points,col='transparent',cex=1.2)

axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)

plot(r~log10(m_abund_reef),data=sp_info,bty='l',xaxt='n',xlab='Mean encounter rate (REEF)',ylab=expression(rho),type='n',ylim=c(-1,1),xlim=c(-3,1.5),cex.lab=1.2)
for(i in 1:nrow(sp_info)){
  lines(c(sp_info$r_l90[i],sp_info$r_u90[i])~rep(log10(sp_info$m_abund_reef[i]),2),lwd=1,col=adjustcolor(col_points[i],alpha.f=0.3))
  lines(rep(sp_info$r[i],2)~c(log10(sp_info$m_abund_reef_l90[i]),log10(sp_info$m_abund_reef_u90[i])),lwd=1,col=adjustcolor(col_points[i],alpha.f=0.3))
}
points(r~log10(m_abund_reef),data=sp_info,pch=21,bg=col_points,col='transparent',cex=1.2)

axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)

plot(r~diff.abund,data=sp_info,bty='l',xlab='Difference (REEF-RVC) in (log10) mean encounter rate',ylab=expression(rho),type='n',ylim=c(-1,1),cex.lab=1.2)
for(i in 1:nrow(sp_info)){
  lines(c(sp_info$r_l90[i],sp_info$r_u90[i])~rep(sp_info$diff.abund[i],2),lwd=1,col=adjustcolor(col_points[i],alpha.f=0.3))
  lines(rep(sp_info$r[i],2)~c(sp_info$diff.abund_l90[i],sp_info$diff.abund_u90[i]),lwd=1,col=adjustcolor(col_points[i],alpha.f=0.3))
}
points(r~diff.abund,data=sp_info,pch=21,bg=col_points,col='transparent',cex=1.2)


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


###sampling sensitivity####

#top 9 species codes:
sp9=c(15,22,60,37,14,66,24,11,44)

#load sensitivity outputs
outs=list.files(here('outputs','sampling effort sensitivity'))
r_sens=list()
for(i in 1:length(outs)){
  r_sens[[i]]=read.csv(here('outputs','sampling effort sensitivity', outs[i]))
  r_sens[[i]]= r_sens[[i]][,-1] #remove row id
}

#plot posteriors for full data and 25% subsamples (20 iterations)
png('FigureS3.png',width=8.5,height=11,units='in',res=300)
par(mfrow=c(3,3))
rho_mode=list()
for(i in 1:length(sp9)){
  plot(c(0,24)~c(-1,1),bty='n',xlab='',ylab='',type='n',yaxt='n',main=sp_info$commonname[sp9[i]])
  if(i==8){
    mtext(side=1,expression(rho),line=3)
  }
  for(q in 1:20){
  d=density(r_sens[[i]][,q],bw=0.02)
  d$y=d$y+(20-q)
  
  polygon(d,col=adjustcolor('blue',alpha.f=0.2),border=NA)
  }
  d2=density(r_list[[sp9[i]]],bw=0.02)
  polygon(d2,col=adjustcolor('darkred',alpha.f=0.5),border=NA)
  rho_mode[[i]]=apply(r_sens[[i]],2,mode)
  
}
dev.off()

png('FigureS4.png',width=8,height=6,units='in',res=300)
plot(c(0.4,1)~c(1,9),ylab='Posterior mode',xlab='Species',bty='l',type='n',xaxt='n')
points(sp_info$mode.r[sp9]~seq(1:9),cex=1.25,bg=adjustcolor('darkred',alpha.f=0.8),pch=21)
for(i in 1:9){
  for(q in 1:20){
    lines(c(sp_info$mode.r[sp9[i]],rho_mode[[i]][q])~c(i,i+0.25))
    
  }
  points(rho_mode[[i]]~rep(i+0.25,20),bg=adjustcolor('blue',alpha.f=0.2),pch=21)
  mtext(sp_info$commonname[sp9[i]],side=1,at=i,cex=0.4)
}
dev.off()

###Variance decomposition####
#For fixed effects (bottom time, depth, etc.) - need variance of effects from data for conditional variance
load(here('data','data.RData'))

files_kl_par_2<- list.files(here('outputs','parameter estimates'))

var_frame_reef_kl<- data.frame(sp_info$commonname,prop.process=NA,prop.obserror=NA,prop.mth=NA,prop.site=NA,prop.hab=NA,prop.strat=NA,prop.diver=NA,prop.d.cluster=NA,prop.m.cluster=NA,prop.btime=NA,prop.depth=NA,prop.vis=NA,prop.current=NA,prop.expert=NA)
var_frame_rvc_kl<- data.frame(sp_info$commonname,prop.process=NA,prop.obserror=NA,prop.mth=NA,prop.psu=NA,prop.hab=NA,prop.strat=NA,prop.depth=NA)
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
  var_full1<- data.frame(sp$sd_q.1.^2,sp$sd_r.1.^2,sp$sd_mth1^2,sp$sd_psu^2,sp$sd_hab1^2,sp$sd_strat1^2,var_beta1)
  var_full2<- data.frame(sp$sd_q.2.^2,sp$sd_r.2.^2,sp$sd_mth1^2,sp$sd_site^2,sp$sd_hab2^2,sp$sd_strat2^2,sp$sd_dv^2,sp$sd_dmy^2,sp$sd_my^2,var_beta2)

  var_frame_rvc_kl[i,2]=median(var_full1[,1])  
  var_frame_rvc_kl[i,3]=median(var_full1[,2])  
  var_frame_rvc_kl[i,4]=median(var_full1[,3])  
  var_frame_rvc_kl[i,5]=median(var_full1[,4])
  var_frame_rvc_kl[i,6]=median(var_full1[,5])  
  var_frame_rvc_kl[i,7]=median(var_full1[,6])
  var_frame_rvc_kl[i,8]=median(var_full1[,7]+var_full1[,8]) #depth + depth2
  
  
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
  var_frame_reef_kl[i,12]=median(var_full2[,11]+var_full2[,12]) #depth + depth2
  var_frame_reef_kl[i,13]=median(var_full2[,13])
  var_frame_reef_kl[i,14]=median(var_full2[,14])
  var_frame_reef_kl[i,15]=median(var_full2[,15])

}
#tables of variance proportions for each factor by species
png('FigureS5.png',width=8.5,height=11,units='in',res=300)
par(mar=c(4,5,2,1),xpd=T)
barplot(height=t(as.matrix(var_frame_rvc_kl[,c(2:8)]/rowSums(var_frame_rvc_kl[,c(2:8)]))),beside=FALSE,horiz=TRUE,
        col=c('#8c2d04',
              '#fdae6b',
              '#fff7bc',
              '#fec44f',
              '#02818a',
              '#016450',
              '#67a9cf'),border=NA)
mtext(side=1,line=2.5,'Proportion of explained variance in RVC survey counts')
mtext(side=2,var_frame_rvc_kl$sp_info.commonname,at=seq(par('usr')[2],par('usr')[4]*0.96,length.out=87),las=2,cex=0.4)
mtext(c('Process var.','Obs var.','Month','PSU','Habitat','Reef Stratum','Depth'),at=c(0.05,0.15,0.25,0.35,0.45,0.55,0.65),cex=0.6,line=-0.5)
points(y=rep(c(par('usr')[4]*1.02),7),x=c(0.05,0.15,0.25,0.35,0.45,0.55,0.65),bg=c(c('#8c2d04',
                                                                    '#fdae6b',
                                                                    '#fff7bc',
                                                                    '#fec44f',
                                                                    '#02818a',
                                                                    '#016450',
                                                                    '#67a9cf')),pch=22)
dev.off()

png('FigureS6.png',width=8.5,height=11,units='in',res=300)
par(mar=c(4,5,2,1),xpd=T)
barplot(height=t(as.matrix(var_frame_reef_kl[,c(2:15)]/rowSums(var_frame_reef_kl[,c(2:15)]))),beside=FALSE,horiz=TRUE,
        col=c('#8c2d04',
              '#fdae6b',
              '#fff7bc',
              '#fec44f',
              '#02818a',
              '#016450',
              '#d9f0a3',
              '#78c679',
              '#238443',
              '#9e9ac8',
              '#67a9cf',
              '#4a1486',
              '#c6dbef',
              '#005a32'),border=NA)
mtext(side=1,line=2.5,'Proportion of explained variance in RVC survey counts')
mtext(side=2,var_frame_rvc_kl$sp_info.commonname,at=seq(par('usr')[2],par('usr')[4]*0.96,length.out=87),las=2,cex=0.4)
mtext(c('Process var.','Obs var.','Month','Site','Habitat','Reef Stratum','Diver ID','Day cluster','Mth. Cluster','Bottom time','Depth','Visibility','Current','Expert'),at=seq(0,1,length.out=14),cex=0.6,line=-0.5)
points(y=rep(c(par('usr')[4]*1.02),14),x=seq(0,1,length.out=14),bg=c('#8c2d04',
                                                                      '#fdae6b',
                                                                      '#fff7bc',
                                                                      '#fec44f',
                                                                      '#02818a',
                                                                      '#016450',
                                                                      '#d9f0a3',
                                                                      '#78c679',
                                                                      '#238443',
                                                                      '#9e9ac8',
                                                                      '#67a9cf',
                                                                      '#4a1486',
                                                                      '#c6dbef',
                                                                      '#005a32'),pch=22)
dev.off()


write.csv(var_frame_rvc_kl,here('outputs','variance decomp','species_variances_rvc_kl.csv'))
write.csv(var_frame_reef_kl,here('outputs','variance decomp','species_variances_reef_kl.csv'))

#Simplified variable categorizations:
#Spatial - Site, Strata
var_frame_reef2=data.frame(procc.var=var_frame_reef_kl[,c('prop.process')],spatial.tot=var_frame_reef_kl[,c('prop.site')],env.tot=rowSums(var_frame_reef_kl[,c('prop.hab','prop.strat','prop.depth')]),observation.tot=rowSums(var_frame_reef_kl[,c('prop.diver','prop.mth','prop.obserror','prop.m.cluster','prop.d.cluster','prop.btime','prop.vis','prop.current','prop.expert')]))

#Simplified variable categorizations:
var_frame_rvc2=data.frame(procc.var=var_frame_reef_kl[,c('prop.process')],spatial.tot=var_frame_rvc_kl[,c('prop.psu')],env.tot=rowSums(var_frame_rvc_kl[,c('prop.hab','prop.strat','prop.depth')]),observation.tot=rowSums(var_frame_rvc_kl[,c('prop.mth','prop.obserror')]))
png('FigureS7.png',width=8.5,height=11,units='in',res=300)
par(mar=c(4,5,2,1),xpd=T)
barplot(height=t(as.matrix(var_frame_rvc2/rowSums(var_frame_rvc2))),beside=FALSE,horiz=TRUE,
        col=c('#8c2d04',
              '#fec44f',
              '#02818a',
              '#fff7bc'),border=NA)
mtext(side=1,line=2.5,'Proportion of explained variance in RVC survey counts')
mtext(side=2,var_frame_rvc_kl$sp_info.commonname,at=seq(par('usr')[2],par('usr')[4]*0.96,length.out=87),las=2,cex=0.4)
mtext(c('Temporal process var.','Spatial','Environmental','Observational'),at=c(0.15,0.3,0.45,0.6),cex=0.6)
points(y=rep(c(par('usr')[4]*1.03),4),x=c(0.15,0.3,0.45,0.6),bg=c('#8c2d04',
                                                                  '#fec44f',
                                                                  '#02818a',
                                                                  '#fff7bc'),pch=22)
dev.off()

png('FigureS8.png',width=8.5,height=11,units='in',res=300)
par(mar=c(4,5,2,1),xpd=T)
barplot(height=t(as.matrix(var_frame_reef2/rowSums(var_frame_reef2))),beside=FALSE,horiz=TRUE,
        col=c('#8c2d04',
              '#fec44f',
              '#02818a',
              '#fff7bc'),border=NA)
mtext(side=1,line=2.5,'Proportion of explained variance in REEF survey counts')
mtext(side=2,var_frame_rvc_kl$sp_info.commonname,at=seq(par('usr')[2],par('usr')[4]*0.96,length.out=87),las=2,cex=0.4)
mtext(c('Temporal process var.','Spatial','Environmental','Observational'),at=c(0.15,0.3,0.45,0.6),cex=0.6)
points(y=rep(c(par('usr')[4]*1.03),4),x=c(0.15,0.3,0.45,0.6),bg=c('#8c2d04',
                                                                  '#fec44f',
                                                                  '#02818a',
                                                                  '#fff7bc'),pch=22)
dev.off()

#tables of category proportions by species
write.csv(var_frame_rvc_kl,here('outputs','variance decomp','species_variancetypes_rvc_kl.csv'))
write.csv(var_frame_reef_kl,here('outputs','variance decomp','species_variancetypes_reef_kl.csv'))


#stacked barplots




