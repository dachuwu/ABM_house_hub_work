# 06/23 cluster plot
require(reshape2);require(ggplot2)
############################################################
################## plotting ################################
setwd("C:/Users/user/Desktop/epidemicABM/biological modeling 105_2/")
source('./ABM_Base_func.R')
source('./ABM_Simulator.R')

result<-RunABM(c( Par.Gamma=(1/6)/24,
                  Par.Beta=10/24,
                  Par.iprob=.05,
                  Par.rel.vuln=1.5,
                  H.total=200,
                  dlen=.5,T.df=1,int.after=0,int.site=5,int.num=5,
                  p.rest=2/3,p.seek=0/3),outName = "test")
#######################################################
GenPlotDat<-function(result,set.var){
  Obs.iState<-result$Obs.iState
  ts.df<-as.data.frame(Obs.iState[,c(2:dim(Obs.iState)[2])])
  ts.df<-cbind(Obs.iState[,1],t(apply(ts.df,1,FUN =function(x) return(c(sum(x==1),sum(x==2),sum(x==3))) )))
  ts.df<-(melt(as.data.frame(ts.df),id.vars = "V1",value.name = "count",variable.name = "state"))
  ts.df$state<-factor(ts.df$state,levels=c("V2","V3","V4"),labels = c("S","I","R"))
  ts.df$set<-set.var
  Obs.iLoc<-result$Obs.iLoc
  loc.ts<-c()
  for (k in c(2:7,0,1)){
    for (i in 1:dim(Obs.iState)[1]){
      sts<-Obs.iState[i,which(Obs.iLoc[i,]==k)]
      loc.ts<-rbind(loc.ts,c(k,Obs.iLoc[i,1],sum(sts==1),sum(sts==2),sum(sts==3)))}
  }
  loc.ts<-as.data.frame(loc.ts)
  colnames(loc.ts)<-c("LOC","time","S","I","R")
  loc.ts.melt<-melt(loc.ts,id.vars=c("LOC","time"))
  loc.ts.melt$LOC<-factor(loc.ts.melt$LOC,levels=c(2:7,0,1),labels = c("U1","U2","U3","C1","C2","C3","Hospital","Household"))
  loc.ts.melt$set<-set.var
  Obs.Inf_Rec<-as.data.frame(result$Obs.Inf_Rec)
  Obs.Inf_Rec$State<-factor(Obs.Inf_Rec$State,levels=c(1:3),labels = c("S","I","R"))
  Obs.Inf_Rec$Loc<-factor(Obs.Inf_Rec$Loc,levels=c(0:7),labels = c("Hospital","Household","U1","U2","U3","C1","C2","C3"))
  Obs.Inf_Rec$set<-set.var
  return(list(ts.df,loc.ts=loc.ts.melt,Obs.Inf_Rec))
}

dataDir<-paste0("SimOutput/0623 set1to6 rep5/",dir(paste0("SimOutput/0623 set1to6 rep5/")))
# s1 7 16 21 23 26
# s2 3 14 24 28 30
# s3 1 8  9 20 29
# s4 4 12 13 17 25
# s5 2 5 11 18 22
# s6 6 10 15 19 27

load(dataDir[1]);#get.info(Result)
PlotDF1<-GenPlotDat(Result,set.var = 1)
load(dataDir[2]);#get.info(Result)
PlotDF2<-GenPlotDat(Result,set.var = 2)
load(dataDir[3]);#get.info(Result)
PlotDF3<-GenPlotDat(Result,set.var = 3)
load(dataDir[4]);#get.info(Result)
PlotDF4<-GenPlotDat(Result,set.var = 4)
load(dataDir[5]);#get.info(Result)
PlotDF5<-GenPlotDat(Result,set.var = 5)
load(dataDir[6]);#get.info(Result)
PlotDF6<-GenPlotDat(Result,set.var = 6)



T.ts.df<-rbind(PlotDF1[[1]],PlotDF2[[1]],PlotDF3[[1]],PlotDF4[[1]],PlotDF5[[1]],PlotDF6[[1]])
T.ts.df$set<-factor(T.ts.df$set,levels = c(1:6))
T.ts.df$V1<-T.ts.df$V1/24
tsplot1<-ggplot(T.ts.df)+theme_bw()+
  geom_path(aes(x=V1,y = count,group=state,color=state))+
  facet_wrap(~set,nrow = 2,ncol = 3)+
  labs(x="time(day)",y="frequency",title="Time series of disease states")+
  annotate("segment",x=1,xend=1,y=-10,yend=700,lty=2,alpha=.8)+
  coord_cartesian(ylim = c(0,650))+
  guides(colour =guide_legend("scenarios"))
ggsave("f5_tsplot1.jpg",plot = tsplot1,scale = .6,width = 18,height = 10)
#
locpattern<-ggplot(PlotDF4[[2]])+
  geom_bar(aes(x = time,y = value,fill=variable),stat="identity",position="stack",width = 1)+
  facet_wrap(~LOC,nrow = 3,ncol = 3,scales = "free_y")+
  theme_classic()
ggsave("f6_locpattern.jpg",plot = locpattern,scale = .6,width = 20,height = 15)

# location of infection events 
T.event.df<-rbind(PlotDF1[[3]],PlotDF2[[3]],PlotDF3[[3]],PlotDF4[[3]],PlotDF5[[3]],PlotDF6[[3]])
T.event.df$set<-factor(T.event.df$set,levels = c(1:6))
case.plot1<-ggplot(data = T.event.df[!is.na(T.event.df$Next)&T.event.df$Next==2,])+
  geom_bar(aes(x=set,fill=Loc),position = "stack")+theme_bw()+
  guides(fill =guide_legend("locale"))+
  labs(x="scenarios")
ggsave("f7_caseplot2.jpg",plot = case.plot1,scale = .6,width = 10,height = 5)

AgeDist1<-ggplot(data = T.event.df[!is.na(T.event.df$Next)&T.event.df$Next==2,])+
  geom_histogram(aes(x=Age,fill=Loc),binwidth = 5,position = "stack")+
  facet_wrap(~set,nrow = 2,ncol = 3)+
  guides(fill =guide_legend("locale"))
ggsave("f8_AgeDist1.jpg",plot = AgeDist1,scale = .6,width = 15,height = 10)


temp<-T.event.df
temp$Age_group<-ifelse(temp$Age<14,"age<=14",ifelse(temp$Age>=65,"age>=65","14<age<65"))
AgeDistTime1<-ggplot(temp[!is.na(temp$Next)&temp$Next==2,])+
  geom_histogram(aes(x = TTE,fill=Age_group),binwidth = 24,position = "fill")+
  facet_wrap(~set,nrow = 2,ncol = 3)
ggsave("AgeDistTime1_non.png",plot = AgeDistTime1,scale = .6,width = 15,height = 10)


### Risk difference

dataDir<-paste0("SimOutput/0623 risk/",dir(paste0("SimOutput/0623 risk/")))
load(dataDir[4])
get.info(Result)
Obs.cumRe<-as.data.frame(Result$Obs.cumRe)
Obs.cumRe$H.sumRe<-rowSums(Obs.cumRe[,c(8:207)],na.rm = T)
dftemp<-Obs.cumRe[,c(1:7,208)]
colnames(dftemp)<-c("time","risk_U1","risk_U2","risk_U3","risk_C1","risk_C2","risk_C3","H.all")
dftemp<-melt(dftemp,id.vars = 'time',variable.name = 'locale',value.name = 'risk')
riskplot1<-ggplot(dftemp)+theme_bw()+
  geom_area(aes(x=time,y=risk,fill=locale))+
  facet_wrap(~locale,nrow = 7)+
  coord_cartesian(ylim=c(0,0.5))
ggsave("riskplot_set1.png",plot = riskplot1,scale = .45,width = 10,height = 20)

load(dataDir[2])
get.info(Result)
Obs.cumRe<-as.data.frame(Result$Obs.cumRe)
Obs.cumRe$H.sumRe<-rowSums(Obs.cumRe[,c(8:207)],na.rm = T)
dftemp<-Obs.cumRe[,c(1:7,208)]
colnames(dftemp)<-c("time","risk_U1","risk_U2","risk_U3","risk_C1","risk_C2","risk_C3","H.all")
dftemp<-melt(dftemp,id.vars = 'time',variable.name = 'locale',value.name = 'risk')
riskplot3<-ggplot(dftemp)+theme_bw()+
  geom_area(aes(x=time,y=risk,fill=locale))+
  facet_wrap(~locale,nrow = 7)+
  coord_cartesian(ylim=c(0,0.5))
ggsave("riskplot_set3.png",plot = riskplot3,scale = .45,width = 10,height = 20)


load(dataDir[5])
get.info(Result)
Obs.cumRe<-as.data.frame(Result$Obs.cumRe)
Obs.cumRe$H.sumRe<-rowSums(Obs.cumRe[,c(8:207)],na.rm = T)
dftemp<-Obs.cumRe[,c(1:7,208)]
colnames(dftemp)<-c("time","risk_U1","risk_U2","risk_U3","risk_C1","risk_C2","risk_C3","H.all")
dftemp<-melt(dftemp,id.vars = 'time',variable.name = 'locale',value.name = 'risk')
riskplot2<-ggplot(dftemp)+theme_bw()+
  geom_area(aes(x=time,y=risk,fill=locale))+
  facet_wrap(~locale,nrow = 7)+
  coord_cartesian(ylim=c(0,0.5))
ggsave("riskplot_set2.png",plot = riskplot2,scale = .45,width = 10,height = 20)


## H.state
load(dataDir[2])
get.info(Result)
Obs.Hstate3<-as.data.frame(Result$Obs.Hstate);
Hsts.df3<-as.data.frame(Obs.Hstate3[,c(2:dim(Obs.Hstate3)[2])])
Hsts.df3<-cbind(Obs.Hstate2[,1],t(apply(Hsts.df3,1,FUN =function(x) return(c(sum(x==1),sum(x==2),sum(x==3))) )))
Hsts.df3<-(melt(as.data.frame(Hsts.df3),id.vars = "V1",value.name = "count",variable.name = "H.state"))
Hsts.df3$H.state<-factor(Hsts.df3$H.state,levels=c("V2","V3","V4"),labels = c("Susceptible","Ongoing","saturated"))
Hsts.df3$set<-3
load(dataDir[5])
get.info(Result)
Obs.Hstate2<-as.data.frame(Result$Obs.Hstate);
Hsts.df2<-as.data.frame(Obs.Hstate2[,c(2:dim(Obs.Hstate2)[2])])
Hsts.df2<-cbind(Obs.Hstate2[,1],t(apply(Hsts.df2,1,FUN =function(x) return(c(sum(x==1),sum(x==2),sum(x==3))) )))
Hsts.df2<-(melt(as.data.frame(Hsts.df2),id.vars = "V1",value.name = "count",variable.name = "H.state"))
Hsts.df2$H.state<-factor(Hsts.df2$H.state,levels=c("V2","V3","V4"),labels = c("Susceptible","Ongoing","saturated"))
Hsts.df2$set<-2
load(dataDir[4])
get.info(Result)
Obs.Hstate1<-as.data.frame(Result$Obs.Hstate);
Hsts.df1<-as.data.frame(Obs.Hstate1[,c(2:dim(Obs.Hstate1)[2])])
Hsts.df1<-cbind(Obs.Hstate1[,1],t(apply(Hsts.df1,1,FUN =function(x) return(c(sum(x==1),sum(x==2),sum(x==3))) )))
Hsts.df1<-(melt(as.data.frame(Hsts.df1),id.vars = "V1",value.name = "count",variable.name = "H.state"))
Hsts.df1$H.state<-factor(Hsts.df1$H.state,levels=c("V2","V3","V4"),labels = c("Susceptible","Ongoing","saturated"))
Hsts.df1$set<-1


Hsts.df<-rbind(Hsts.df1,Hsts.df2,Hsts.df3)
Hsts.df$set<-factor(Hsts.df$set,levels = c(1:3))
Hsts.df$V1<-Hsts.df$V1/24

h.stsplot<-ggplot(Hsts.df)+theme_bw()+
  geom_path(aes(x=V1,y = count,group=H.state,color=H.state))+
  facet_wrap(~set,nrow = 3,ncol = 1)+
  labs(x="time(day)",y="frequency",title="Time series of household states")

ggsave("h.stsplot.png",plot = h.stsplot,scale = .6,width = 10,height = 15)


Ags<-create.pop(paras = def.set)
View(Result$Vinces)
