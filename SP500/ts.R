max_carmdd<-function(data,para){
    final = 500
    tpf = seq(0.1,2,by=.1)
    slf = seq(0.1,2,by=.1)
    carmdd = matrix(0,nrow=length(tpf),ncol=length(slf))
    for(i in 1:length(tpf))
    for(j in 1:length(slf)){
        sig = get_sig(para=para,from=1,to=final)
        p = trade_stats1(signal=sig,data=data,sl_fac = slf[j],tp_fac=tpf[i])
        eq = p$eq
        eq = eq[2:length(eq)]/eq[1:length(eq)-1]-1
        t = xts(eq,index(data)[2:500])
        if(!is.na(coredata(Return.annualized(t)/maxDrawdown(t))))
        carmdd[i,j] = coredata(Return.annualized(t)/maxDrawdown(t))
    }
    
    
    return(carmdd)
}

max_carmdd_vol<-function(data,para){
    final = 500
    lb_vol = seq(1,20,by=1)
    roll = seq(1,15,by=1)
    carmdd = matrix(0,nrow=length(lb_vol),ncol=length(roll))
    for(i in 1:length(lb_vol))
    for(j in 1:length(roll)){
        sig = get_sig_vol(para=para,data=data,vol_lb=lb_vol[i],rolling_vol=roll[j],from=1,to=final)
        p = trade_stats1(signal=sig,data=data)
        eq = p$eq
        eq = eq[2:length(eq)]/eq[1:length(eq)-1]-1
        t = xts(eq,index(data)[2:500])
        if(!is.na(coredata(Return.annualized(t)/maxDrawdown(t))))
        carmdd[i,j] = Return.annualized(t)/maxDrawdown(t)
    }
    
    
    
    return(carmdd)
}
get_sig<-function(para,threshold=0,from=1,to=NA){
    mu = c(rep(0,29),para$mu)
    n = length(mu)
    dir_mu = rep(0,n)
    dir_mu[mu>threshold] = 1
    dir_mu[mu<(-threshold)] = -1
    bb = BBands(Cl(GSPC),maType='EMA',n=10,sd=2)
    bb_sig = rep(0,length(Cl(GSPC)))
    bb_sig[bb$mavg<Cl(GSPC)] = -1
    bb_sig[bb$mavg>Cl(GSPC)] = 1
    dir = rep(0,length(Cl(GSPC)))
    dir[bb_sig+dir_mu==2] = 1
    dir[bb_sig+dir_mu==-2] = -1
    dir = c(0,Lag(dir)[-1])
    if(!is.na(to))
        dir = dir[from:to]
    return(dir)
}

get_sig_vol<-function(para,threshold=0,data,vol_lb=5,vol_filter=0.5,rolling_vol=1,from=1,to=NA){
    mu = c(rep(0,29),para$mu)
    n = length(mu)
    sig = rep(0,n)
    ret = ((Cl(data) - Op(data))/Op(data)*100)
    vol = rep(0,n)
    vol[ret>vol_filter] = 1
    vol[ret<(-vol_filter)]=-1
    sum_vol=filter(vol,rep(1,vol_lb), sides=1)
    sig_vol = rep(0,n)
    sig_vol[sum_vol>=rolling_vol] = -1
    sig_vol[sum_vol<=(-rolling_vol)] = 1
    dir_mu = rep(0,n)
    dir_mu[mu>threshold] = 1
    dir_mu[mu<(-threshold)] = -1
    sig[(dir_mu+sig_vol)==2] = 1
    sig[(dir_mu+sig_vol)==-2] = -1
    sig = c(0,Lag(sig)[-1])
    if(!is.na(to))
    sig = sig[from:to]
    return(sig)
}

get_trans<-function(signal,data,cost=0.005){
      trans=rep(1,length(signal))
      trans[signal==0] = 0
      for(i in 2:length(signal)){
            if(signal[i] !=0)
                  if(coredata(signal[i])==coredata(signal[i-1]))
                        trans[i] = 0
      }
      return(cost*trans)
}
trade_stats_b<-function(signal,data,starting=10000,transcost=0.005,delta=0.02,point = 5,sl_fac=1.0,tp_fac=1.0){
    curr_pos = 0
    lots = 1
    curr_eq = starting
    dir = c()
    before=c()
    after=c()
    profit=c()
    lot=c()
    entry=c()
    sl=c()
    length=c()
    pc_profit=c()
    tp=c()
    sl_pt=c()
    eq=c()
    curr = 1
    flag = 1
    for(i in 1:length(signal)){
          if(curr_pos==0){
                if(signal[i]==-1){
                      curr_pos = 1
                      lot[curr] = 1
                      curr_eq = curr_eq - coredata(transcost*Op(data)[i]*lot[curr])
                      dir[curr] = 1
                      before[curr] = curr_eq
                      entry[curr] = Op(data)[i]
                      length[curr] = 1
                }

          }
          else{
                length[curr] = length[curr] + 1
                if(signal[i] !=1){
                      sl[curr] = Op(data)[i]
                      profit[curr] = (sl[curr] - entry[curr])*dir[curr]*lot[curr]
                      after[curr] = profit[curr] + before[curr]
                      pc_profit[curr] = (sl[curr] - entry[curr])/entry[curr]
                      eq[curr] = after[curr]
                  
                      curr = curr+1
                      curr_pos=0
                }
          }
    }
    t_s = list(dir=dir,before=before,after=after,profit=profit,lot=lot,entry=entry,sl=sl,length=length,pc_profit=pc_profit,tp=tp,sl_pt=sl_pt,eq=eq)

    return(t_s)
}