library(quantmod)
library(PerformanceAnalytics)
library(plotly)
getSymbols('^GSPC',src='yahoo',from='2007-1-3',to='2016-3-4')

max_carmdd<-function(data,para){
      final = 500
      lb_vol = seq(2,25,by=1)
      roll = seq(2,25,by=1)
      carmdd = matrix(0,nrow=length(lb_vol),ncol=length(roll))
      for(i in 1:length(lb_vol))
            for(j in 1:length(roll)){
                  sig = get_sig_vol(para=para,data=data,vol_lb=lb_vol[i],rolling_vol=roll[j],from=1,to=final)
                  p = trade_stats(signal=sig,data=data)
                  eq = p$eq
                  eq = eq[2:length(eq)]/eq[1:length(eq)-1]-1
                  t = xts(eq,index(GSPC)[2:500])
                  if(!is.na(coredata(Return.annualized(t)/maxDrawdown(t))))
                        carmdd[i,j] = coredata(Return.annualized(t)/maxDrawdown(t)) 
            }
      
      
      return(carmdd)
}
max_carmdd_trail<-function(data,para){
      final = 500
      trail = seq(1,10,by=1)
      sl_fac = seq(0.1,2,by=0.1)
      carmdd = matrix(0,nrow=length(trail),ncol=length(sl_fac))
      for(i in 1:length(trail))
            for(j in 1:length(sl_fac)){
                  sig = get_sig_vol(para=para,data=data,vol_lb=7,rolling_vol=2,from=1,to=final)
                  p = trade_stats(signal=sig,data=data,trail=trail[i],sl_fac=sl_fac[j])
                  eq = p$eq
                  eq = eq[2:length(eq)]/eq[1:length(eq)-1]-1
                  t = xts(eq,index(GSPC)[2:500])
                  if(maxDrawdown(t)<0.1)
                        if(!is.na(coredata(Return.annualized(t)/maxDrawdown(t))))
                              carmdd[i,j] = coredata(Return.annualized(t)/maxDrawdown(t)) 
            }
      
      
      return(carmdd)     
}

get_SP500_data<-function(lookback){
      getSymbols("^GSPC",src='yahoo')
      n = length(Cl(GSPC))
      d = as.data.frame(log(coredata(Cl(GSPC))))
      for(i in lookback:n){
            write.csv(d$GSPC.Close[(i-lookback+1):i],file=sprintf('data/data%d.csv',(i-lookback+1)))
      }    
}

get_graphs<-function(n){
      mu = c()
      gam = c()
      sigma_j = c()
      sigma_v = c()
      rho = c()
      k = c()
      v = c()
      for(i in 1:n){
            dat = read.csv(sprintf('output/output%d.csv',i))
            mu[i] = sum(dat$mu*dat$norm_w)
            gam[i] = sum(dat$gam*dat$norm_w)
            sigma_j[i] = sum(dat$sigma_j*dat$norm_w)
            rho[i] = sum(dat$rho*dat$norm_w)
            k[i] = sum(dat$k*dat$norm_w)
            v[i] = sum(dat$v*dat$norm_w)
            sigma_v[i] = sum(dat$sigma_v*dat$norm_w)
      }
      return(list(mu=mu,gam=gam,sigma_j=sigma_j,sigma_v = sigma_v,rho=rho,k=k,v=v))
}

get_sig<-function(para,threshold=0){
      n = length(para$mu)
      sig = rep(0,n)
      add = rep(0,30)
      sig[para$mu>threshold] = 1
      sig[para$mu<(-threshold)] = -1
      sig = c(add,Lag(sig)[-1])
      return(sig)
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

trade_stats<-function(signal,data,starting=10000,transcost=0.005,trail=3,delta=0.02,point = 5,sl_fac=1.0){
      curr_pos = 0
      lots = 1
      curr_eq = starting
      eq_curve = starting
      eq_curve_low = starting
      t_s = list(dir=c(),before=c(),after=c(),profit=c(),lot=c(),entry=c(),sl=c(),length=c(),pc_profit=c())
      curr = 1
      flag = 1
      for(i in 1:length(signal)){
      	if(i>1){
      		eq_curve[i] = eq_curve[i-1]
      		eq_curve_low[i] = eq_curve_low[i-1]
      	}
            if(curr_pos !=0){
                  t_s$length[curr] = t_s$length[curr]+1
                  if(t_s$dir[curr]==1){
                        if(b_trailing>coredata(Lo(data)[i])){
                              t_s$sl[curr] = b_trailing
                              t_s$profit[curr] = (t_s$sl[curr] - t_s$entry[curr])*t_s$lot[curr]*t_s$dir[curr]*point
                              t_s$pc_profit[curr] = t_s$profit[curr]/t_s$before[curr]
                              eq_curve_low[i] = curr_eq + (b_trailing - t_s$entry[curr])*t_s$lot[curr]*t_s$dir[curr]*point
                              curr_eq = curr_eq+t_s$profit[curr]
                              t_s$after[curr] = curr_eq
                              curr_pos = 0
                              eq_curve[i] = eq_curve[i] + (b_trailing - coredata(Op(data)[i-1]))*t_s$lot[curr]*t_s$dir[curr]*point
                              
                              curr = curr+1   	
                              
                        }
                        else{
                              b_trailing = max(b_trailing,min(coredata(Lo(data)[i-trail+1:i])))
                              eq_curve[i] = eq_curve[i] + coredata(momentum(Op(data))[i])*t_s$lot[curr]*t_s$dir[curr]*point
                              eq_curve_low[i] = curr_eq+(b_trailing-t_s$entry[curr])*t_s$lot[curr]*t_s$dir[curr]*point
                              

                        }
                  }
                  else{
                        if(t_s$dir[curr]==-1){
                              if(s_trailing<coredata(Hi(data)[i])){
                                    t_s$sl[curr] = s_trailing
                                    t_s$profit[curr] = (t_s$sl[curr] - t_s$entry[curr])*t_s$lot[curr]*t_s$dir[curr]*point
                                    t_s$pc_profit[curr] = t_s$profit[curr]/t_s$before[curr]
                                    eq_curve_low[i] = curr_eq+(s_trailing-t_s$entry[curr])*t_s$lot[curr]*t_s$dir[curr]*point
                                    curr_eq = curr_eq+t_s$profit[curr]
                                    t_s$after[curr] = curr_eq
                                    curr_pos = 0
                                    eq_curve[i] = eq_curve[i]+(s_trailing-coredata(Op(data)[i-1]))*t_s$lot[curr]*t_s$dir[curr]*point
                                    
                                    
                                    curr = curr+1
                              }
                              else{
                                    s_trailing = min(s_trailing,max(coredata(Hi(data)[i-trail+1:i])))
                                    eq_curve[i] = eq_curve[i] + coredata(momentum(Op(data))[i])*t_s$lot[curr]*t_s$dir[curr]*point
						eq_curve_low[i]=curr_eq+(s_trailing-t_s$entry[curr])*t_s$lot[curr]*t_s$dir[curr]*point
									

                              		
                              }
                        }
                  }

            }
            else{
                  if(signal[i]==1){
                        b_trailing = min(coredata(Lo(data))[i-trail+1],coredata(Op(data)[i])-sl_fac*(coredata(Hi(data))[i-1]-coredata(Lo(data))[i-1]))
                        t_s$entry[curr] = coredata(data)[i,1]
                        t_s$lot[curr] = max(1,round(curr_eq*delta/(point*(t_s$entry[curr]-b_trailing))))
                        t_s$dir[curr] = 1
                        curr_eq = curr_eq - transcost*t_s$lot[curr]*t_s$entry[curr]
                        t_s$before[curr] = curr_eq
                        
                        t_s$length[curr] = 1
                        curr_pos = 1
                        eq_curve[i] = eq_curve[i]- transcost*t_s$lot[curr]*t_s$entry[curr]
                        eq_curve_low[i] = curr_eq+(b_trailing-t_s$entry[curr])*t_s$lot[curr]*t_s$dir[curr]*point
                        
                        
                  }
                  else{
                        if(signal[i]==-1){
                              s_trailing = max(coredata(Hi(data))[i-trail+1],sl_fac*(coredata(Hi(data))[i-1]-coredata(Lo(data))[i-1])+coredata(Op(data))[i])
                              t_s$entry[curr] = coredata(data)[i,1]
                              t_s$lot[curr] = max(1,round(curr_eq*delta/(point*(s_trailing-t_s$entry[curr]))))
                              t_s$dir[curr] = -1
                              curr_eq = curr_eq - transcost*t_s$lot[curr]*t_s$entry[curr]
                              t_s$before[curr] = curr_eq
                              t_s$length[curr] = 1
                              curr_pos = -1
                              eq_curve[i] = eq_curve[i]- transcost*t_s$lot[curr]*t_s$entry[curr]
					eq_curve_low[i] = curr_eq+(s_trailing-t_s$entry[curr])*t_s$lot[curr]*t_s$dir[curr]*point
                         	  
                        }
                  }
            }
      }
      i = length(signal)
      if(curr_pos !=0){
          t_s$sl[curr] = coredata(Cl(data)[i])
          t_s$profit[curr] = (t_s$sl[curr] - t_s$entry[curr])*t_s$lot[curr]*t_s$dir[curr]*point
          t_s$pc_profit[curr] = t_s$profit[curr]/t_s$before[curr]
          eq_curve_low[i] = curr_eq + (t_s$sl[curr] - t_s$entry[curr])*t_s$lot[curr]*t_s$dir[curr]*point
          curr_eq = curr_eq+t_s$profit[curr]
          t_s$after[curr] = curr_eq
          curr_pos = 0
          eq_curve[i] = eq_curve[i] + (b_trailing - coredata(Op(data)[i-1]))*t_s$lot[curr]*t_s$dir[curr]*point
      }

      t_s$low = eq_curve_low
      t_s$eq = eq_curve
      return(t_s)
}