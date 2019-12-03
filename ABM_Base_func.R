### base functions
States<-c(1:3) 
draw.age<-function()sample(c(1:85),size = 1,prob = rep(c(14,74,12),times=c(14,50,21)))# 14-74-12
draw.route<-function(age){
  if(age<=14) { sample(c(1:9),size = 1,prob = c(8,1,1,8,1,1,8,1,1))
  } else if (age<=65){ sample(c(1:9),size = 1,prob = c(1,8,1,1,8,1,1,8,1))
  } else  sample(c(1:9),size = 1,prob = c(1,1,8,1,1,8,1,1,8) )
}
Route_mtx<-matrix(c(rep(1,9),rep(2:4,each=3),rep(5:7,3),rep(2:4,each=3)),9,4)
# route1: U1-C1  1 2 5 2
# route2: U1-C2  1 2 6 2
# route3: U1-C3  1 2 7 2
# route4: U2-C1  1 3 5 3
# route5: U2-C2  1 3 6 3 
# route6: U2-C3  1 3 7 3 
# route7: U3-C1  1 4 5 4 
# route8: U3-C2  1 4 6 4 
# route9: U3-C3  1 4 7 4 
# route10: hospitalization 0 1 
new_agent <- function(id, st=States[1], age_rdist=draw.age, route_rdist=draw.route, loc=1,paras){
  with(as.list(paras),{
  age<-age_rdist()
  return (c(ID=id,H.ID=NA,Age=age, Route=route_rdist(age), 
            log.Loc=1, Loc=loc, nex.Loc=NA,TTM=Inf,Rate.move=0,
            State=st, Next=NA, TTE=Inf, Risk=0,Seek=0,
            P.inf=if(age>14&age<65) Par.iprob*(1-(Par.rel.vuln*.26))/.74 else Par.iprob*Par.rel.vuln) )
  })
}
create.pop<-function(H.size.dist=c(1629970, 1633310, 1551340, 1412485, 647519, 539747),paras){
  with(as.list(paras),{
      H.size<-sample(c(1:6),prob=H.size.dist,H.total,replace = T)
  H.id<-rep(c(1:H.total),times=H.size)
  n.pop<-length(H.id)
  Vinces<-c() 
  for(i in 1:n.pop){Vinces<-rbind(Vinces,new_agent(id = i,paras = paras))}
  Vinces[,'H.ID']<-sample(x = H.id,size = length(H.id),replace = F)
  return(Vinces)
  })
}

Find.nxt.event <- function(Agent,Time=0,paras){
  with(as.list(paras),{
    if(Agent['State']==States[1]){
    if(Agent['Risk']>0){
      Agent['Next']<-States[2]
      Agent['TTE']<-rexp(1,Agent['Risk']) + Time
    } 
  } else if(Agent['State']==States[2]){
    Agent['Next']<-States[3]
    Agent['TTE']<-rexp(1,Par.Gamma) + Time
  }
  return(Agent)
  })
}

Find.nxt.Loc<-function(Agent,Time=0,route_mtx=Route_mtx){
  i <- Agent['Route'] 
  j <- Agent['log.Loc'] %% 4
  Agent['nex.Loc']<-route_mtx[i,(j+1)]
  t.int<- Time %% 24 
  k<-.1
  l<-(12-k)/11
  
  if(Agent['State']!=2){
    if( (t.int>=14&t.int<15)|(t.int>=23&t.int<24) ){
      if(Agent['Loc']==1) {
        Agent['TTM']<-rexp(1,1/(14*k))+Time
      } else if (Agent['Loc']%in%c(0,2:4)) {
        Agent['TTM']<-rexp(1,1/1)+Time
      } else if (Agent['Loc']%in%c(5:7)) {
        Agent['TTM']<-rexp(1,1/(8*k))+Time} 
    } else {
      if(Agent['Loc']==1) {
        Agent['TTM']<-rexp(1,1/(14*l))+Time
      } else if (Agent['Loc']%in%c(0,2:4)) {
        Agent['TTM']<-rexp(1,1/1)+Time
      } else if (Agent['Loc']%in%c(5:7)) {
        Agent['TTM']<-rexp(1,1/(8*l))+Time}
    }
  } else if (Agent['State']==2){
    if(Agent['Seek']==1){
      if(Agent['Loc']==1) {
        Agent['TTM']<-rexp(1,1/21)+Time
      } else if (Agent['Loc']%in%c(2:4)) {
        Agent['TTM']<-rexp(1,1/1)+Time
      } else if (Agent['Loc']%in%c(5:7)) {
        Agent['TTM']<-rexp(1,1/1)+Time}
      
    } else if (Agent['Seek']==2){
      Agent['nex.Loc']<-0
      Agent['TTM'] <- if(Agent['Loc']!=0) rexp(1,1/24)+Time else Inf
      
    } else {
      if( (t.int>=14&t.int<15)|(t.int>=23&t.int<24) ){
        if(Agent['Loc']==1) {
          Agent['TTM']<-rexp(1,1/(14*k))+Time
        } else if (Agent['Loc']%in%c(2:4)) {
          Agent['TTM']<-rexp(1,1/1)+Time
        } else if (Agent['Loc']%in%c(5:7)) {
          Agent['TTM']<-rexp(1,1/(8*k))+Time}
      } else {
        if(Agent['Loc']==1) {
          Agent['TTM']<-rexp(1,1/(14*l))+Time
        } else if (Agent['Loc']%in%c(2:4)) {
          Agent['TTM']<-rexp(1,1/1)+Time
        } else if (Agent['Loc']%in%c(5:7)) {
          Agent['TTM']<-rexp(1,1/(8*l))+Time}
      }
    }
  }
  return(Agent)
}

transit.St <-function(Ags,pos){
  Ags[pos,'State']<-Ags[pos,'Next']
  Ags[pos,'TTE']<-Inf
  return(Ags)
}
transit.Loc<-function(Ags,pos){
  Ags[pos,'Loc']<-Ags[pos,'nex.Loc']
  Ags[pos,'log.Loc']<-Ags[pos,'log.Loc']+1
  Ags[pos,'TTM']<-Inf
  return(Ags)
}
cal.risk<-function(Ags,paras){
  with(as.list(paras),{
    for(j in unique(Ags[which(Ags[,'Loc']==1),'H.ID'])){
    Ags[which(Ags[,'Loc']==1&Ags[,'H.ID']==j),'Risk']<-mean(Ags[which(Ags[,'Loc']==1&Ags[,'H.ID']==j),'State']==States[2])
  }
  for(i in 2:7){
    Ags[which(Ags[,"Loc"]==i),'Risk']<-mean(Ags[which(Ags[,"Loc"]==i),'State']==States[2])
  }
  Ags[,'Risk']<- Ags[,'Risk']*Par.Beta*Ags[,'P.inf']
  return(Ags)
  })
}

## Observe function
observe.State <- function(Ags, time){
  #obs <- c(table(factor(Ags[, 'State'], levels=0:2)))
  return (c(Time=time, Ags[,'State']))
}
observe.Loc<-function(Ags, time){
  return (c(Time=time, Ags[,'Loc']))
}
observe.event<-function(df,Ags,e){
  return(rbind(df,Ags[e,]))
}
observe.Re<-function(Ags,time,paras){
  with(as.list(paras),{
    H.Re<-c()
    for (k in c(1:H.total)){
      H.Re<-c(H.Re,
              Par.Beta*
                mean(Ags[which(Ags[,"Loc"]==1&Ags[,"H.ID"]==k&Ags[,'State']==1),'P.inf'])*
                mean(Ags[which(Ags[,"Loc"]==1&Ags[,"H.ID"]==k),'State']==1)*length(which(Ags[which(Ags[,"Loc"]==1&Ags[,"H.ID"]==k),'State']==2)) )
    }
    L.Re<-c()
    for(i in c(2:7)){
      L.Re<-c(L.Re,
              Par.Beta*
                mean(Ags[which(Ags[,"Loc"]==i),'P.inf'])*
                mean(Ags[which(Ags[,"Loc"]==i),'State']==1)*
                length(which(Ags[which(Ags[,"Loc"]==i),'State']==2)) )
    }
    seq.Re<-c(Time=time,L.Re,H.Re)
    #names(sumRe)<-c("time",paste0("L.Re",c(2:7)),paste0("H.Re",c(1:H.total)))
    return(seq.Re)
  })
}

calc.Hstate<-function(x){
  if(sum(x!=0)==0){
    return(0)
  }else if(x[1]==0&x[2]==0&x[3]!=0){
    return(3)
  }else if(x[2]!=0){
    return(2)
  }else if(x[1]!=0){
    return(1)
  }
}

observe.Hstate<-function(Ags,time){
  htab<-c()
  for(i in 1:max(Ags[,'H.ID'])){
    htab<-rbind(htab,
                c(S=length(which(Ags[,'H.ID']==i&Ags[,'State']==1)),
                  I=length(which(Ags[,'H.ID']==i&Ags[,'State']==2)),
                  R=length(which(Ags[,'H.ID']==i&Ags[,'State']==3)))
    )}
  return(c(Time=time,apply(htab,MARGIN = 1,calc.Hstate)))
}



get.info_o<-function(SimOutput=list()){
  with(SimOutput,{
    message_run<-(paste0("Introduction: ",int.num," infectives at site ",int.site," after ",int.after/24," day elapsed\n",
                         "Adaptation: [none] = ",round(c(1-p.seek-p.rest)*100,digits=1),"%\n",
                         "Adaptation: [rest] = ",round(c(p.rest)*100,digits=1),"%\n",
                         "Adaptation: [seek] = ",round(c(p.seek)*100,digits=1),"%\n",
                         "Par.effective contact rate = ",Par.Beta*24," per day\n",
                         "Par.mean duration of infectiousness = ",(Par.Gamma*24)^-1," day\n",
                         "H.total: ",H.total,"\nno.individuals: ",dim(Vinces)[1],
                         "\nSimulation length: ",dlen," days","\nRecord diff: ",T.df," hr\nAge-assortitive movement."))
    cat(message_run)
  })
}
get.info<-function(SimOutput=list()){
  with(SimOutput,{
    message_run<-(paste0("Introduction: ",int.num," infectives at site ",int.site," after ",int.after/24," day elapsed\n",
                         "Adaptation: [none] = ",round(c(1-p.seek-p.rest)*100,digits=1),"%\n",
                         "Adaptation: [rest] = ",round(c(p.rest)*100,digits=1),"%\n",
                         "Adaptation: [seek] = ",round(c(p.seek)*100,digits=1),"%\n",
                         "Par.effective contact rate = ",Par.Beta*24," per day\n",
                         "Par.probability of infection = ",round(Par.iprob*100,digits=1),"%\n",
                         "Par.relative susceptibility = ",Par.rel.vuln," \n",
                         "Par.mean duration of infectiousness = ",(Par.Gamma*24)^-1," day\n",
                         "H.total: ",H.total,"\nno.individuals: ",dim(Vinces)[1],
                         "\nSimulation length: ",dlen," days","\nRecord diff: ",T.df," hr\nAge-assortitive movement."))
    cat(message_run)
  })
}
######################################################################
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