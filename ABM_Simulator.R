#source('~/Documents/R project/biological modeling 105_2/ABM_Base_func.R')
def.set<-c( Par.Gamma=(1/6)/24,
            Par.Beta=10/24,
            Par.iprob=.05,
            Par.rel.vuln=1.5,
            H.total=200,
            dlen=30,
            T.df=1,
            int.after=0,
            int.site=5,
            int.num=5,
            p.rest=1/3,p.seek=1/3)

with(as.list(def.set),{Par.Beta*24})
RunABM<-function(Settings=def.set,outName=""){
  with(as.list(Settings),{
  T.start<-0
  T.end<-24*dlen
  T.now<-T.start
  T.rec<-T.start
  ini<-0
  Vinces<-create.pop(paras = Settings)
  Vinces<-cal.risk(Vinces,paras = Settings)
  Vinces<-t(apply(Vinces, 1, Find.nxt.event,Time = T.now,paras=Settings))
  Vinces<-t(apply(Vinces, 1, Find.nxt.Loc,Time = T.now))
  # observe functions
  Obs.iState<-observe.State(Vinces,T.now)
  Obs.iLoc<-observe.Loc(Vinces,T.now)
  Obs.cumRe<-observe.Re(Vinces,T.now,paras = Settings)
  Obs.Hstate<-observe.Hstate(Vinces,T.now)
  Obs.Inf_Rec<-c()
  ############ Simulator begin #########################################
  pb <- txtProgressBar(min = 0, max = T.end, initial = 0,style = 3 )
  
  
  while(T.now<T.end){
    while(T.rec+T.df < T.now){
      T.rec <- T.rec+T.df
      Obs.iState<-rbind(Obs.iState,observe.State(Vinces,time = T.rec))
      Obs.iLoc<-rbind(Obs.iLoc,observe.Loc(Vinces,T.rec))
      Obs.cumRe<-rbind(Obs.cumRe,observe.Re(Vinces,T.rec,paras = Settings))
      Obs.Hstate<-rbind(Obs.Hstate,observe.Hstate(Vinces,T.rec))
    }
    if (min(Vinces[,'TTM'])<min(Vinces[,'TTE'])){
      nxt.m<-which.min(Vinces[,'TTM'])
      if(Vinces[nxt.m,'TTM']>=T.end) break
      T.now<-Vinces[nxt.m,'TTM']
      Vinces<-transit.Loc(Vinces,nxt.m)
      if(Vinces[nxt.m,'Loc']==0) cat(paste0("\n",T.now," Event : Hospitalization\n"))
        if(T.now>int.after & ini<int.num & Vinces[nxt.m,'Loc']==int.site){
          Vinces[nxt.m,'State']<-2
          Vinces[nxt.m,'Seek']<-sample(c(0,1,2),size = 1,prob = c(1-p.seek-p.rest,p.rest,p.seek))
          cat(paste0("\n",T.now," Event : 1 infective introduced\n"))
          Obs.Inf_Rec<-observe.event(Obs.Inf_Rec,Vinces,nxt.m)
          ini<-ini+1
        }
    } else if (min(Vinces[,'TTM'])>min(Vinces[,'TTE'])){
      nxt.e<-which.min(Vinces[,'TTE'])
      if(Vinces[nxt.e,'TTE']>=T.end) break
      T.now<-Vinces[nxt.e,'TTE']
      Obs.Inf_Rec<-observe.event(Obs.Inf_Rec,Vinces,nxt.e)
      Vinces<-transit.St(Vinces,nxt.e)
      if(Vinces[nxt.e,'State']==2) {
        Vinces[nxt.e,'Seek']<-sample(c(0,1,2),size = 1,prob = c(1-p.seek-p.rest,p.rest,p.seek))
        cat(paste0("\n",T.now," Event : infection\n"))
      }
    } else {
      nxt.m<-which.min(Vinces[,'TTM'])
      if(Vinces[nxt.m,'TTM']>=T.end) break
      T.now<-Vinces[nxt.m,'TTM']
      Vinces<-transit.Loc(Vinces,nxt.m)
        if(T.now>int.after & ini<int.num & Vinces[nxt.m,'Loc']==int.site){
          Vinces[nxt.m,'State']<-2
          cat(paste0("\n",T.now," Event : 1 infective introduced\n"))
          Obs.Inf_Rec<-observe.event(Obs.Inf_Rec,Vinces,nxt.m)
          ini<-ini+1
        }
      nxt.e<-which.min(Vinces[,'TTE'])
      Obs.Inf_Rec<-observe.event(Obs.Inf_Rec,Vinces,nxt.e)
      Vinces<-transit.St(Vinces,nxt.e)
      if(Vinces[nxt.e,'State']==2) {
        Vinces[nxt.e,'Seek']<-sample(c(0,1,2),size = 1,prob = c(1-p.seek-p.rest,p.rest,p.seek))
        cat(paste0("\n",T.now," Event : infection\n"))
      }
    }
    Vinces<-cal.risk(Vinces,paras=Settings)
    Vinces<-t(apply(Vinces, 1, Find.nxt.event,Time = T.now,paras=Settings))
    Vinces<-t(apply(Vinces, 1, Find.nxt.Loc,Time = T.now))
    setTxtProgressBar(pb, T.now)
  }
  
  while(T.rec+T.df < T.end){
    T.rec <- T.rec+T.df
    Obs.iState<-rbind(Obs.iState,observe.State(Vinces,time = T.rec))
    Obs.iLoc<-rbind(Obs.iLoc,observe.Loc(Vinces,T.rec))
    Obs.cumRe<-rbind(Obs.cumRe,observe.Re(Vinces,T.rec,paras = Settings))
    Obs.Hstate<-rbind(Obs.Hstate,observe.Hstate(Vinces,T.rec))
  }
  
  if(T.rec+T.df >= T.end){
    Obs.iState<-rbind(Obs.iState,observe.State(Vinces,time = T.end))
    Obs.iLoc<-rbind(Obs.iLoc,observe.Loc(Vinces,T.end))
    Obs.cumRe<-rbind(Obs.cumRe,observe.Re(Vinces,T.end,paras = Settings))
    Obs.Hstate<-rbind(Obs.Hstate,observe.Hstate(Vinces,T.rec))
  }
  ############ Simulator end ###########################################
  Result<-list(Obs.iState=Obs.iState,Obs.iLoc=Obs.iLoc,Obs.Inf_Rec=Obs.Inf_Rec,
               Vinces=Vinces,Obs.cumRe=Obs.cumRe,Obs.Hstate=Obs.Hstate,
               dlen=dlen,H.total=H.total,int.after=int.after,int.num=int.num,int.site=int.site,
               p.rest=p.rest,p.seek=p.seek,Par.Beta=Par.Beta,Par.Gamma=Par.Gamma,
               Par.iprob=Par.iprob,Par.rel.vuln=Par.rel.vuln,
               T.df=T.df)
  save(Result,file = paste0("SimOutput/ABMoutput_",outName,"_",format(Sys.time(),"%Y%m%d_%H_%M_%S"),".RDat"))
  return(Result)
  })
  
}
