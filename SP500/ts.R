library(quantmod)
library(PerformanceAnalytics)
library(plotly)
getSymbols('^GSPC',src='yahoo',from='2007-1-3',to='2016-3-4')

get_graphs<-function(n=2280){
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

get_trans<-function(signal,data,cost=0.001){
      trans=rep(1,length(signal))
      trans[signal==0] = 0
      for(i in 2:length(signal)){
            if(signal[i] !=0)
                  if(coredata(signal[i])==coredata(signal[i-1]))
                        trans[i] = 0
      }
      return(cost*trans)
}

weighted_portfolio<-function(s1,s2,s3,s4,from=1,to=100){
      w1 = seq(0,1,by=0.1)
      w2 = seq(0,1,by=0.1)
      w3 = seq(0,1,by=0.1)
      output = c()
      for(i in w1)
            for(j in w2)
                  for(k in w3){
                        w = 1-i-j-k
                        if(w>0){
                              comb = i*diff(s1)+j*diff(s2)+k*diff(s3)+w*diff(s4)
                              comb[1] = 0
                              stat1 = Return.annualized(comb[from:to])
                              stat = stat1/maxDrawdown(comb[from:to])
                              output = rbind(output,c(i,j,k,w,stat,stat1)) 
                        }

                  }
      return(output)
}
rolling_opt<-function(s1,s2,s3,s4,n=50){
      time = seq(100,2309-n,by=n)
      k = round(100/n)
      p = (0.6*diff(s1) + 0.1*diff(s2) +0*diff(s3)+0.3*diff(s4))[101:200]
      weights = c(0.6,0.1,0,0.3)
      for(i in (k+1):(length(time))){
            opt_from = time[i-k]
            opt_to = time[i]
            w = weighted_portfolio(s1,s2,s3,s4,from=opt_from,to=opt_to)
            w = w[order(w[,5],w[,6],decreasing=T),]
            if(is.na(w[1,5])){
                  w1 = 0.25
                  w2 = 0.25
                  w3 = 0.25
                  w4 = 0.25
            }
            else{
                  w1 = w[1,1]
                  w2 = w[1,2]
                  w3 = w[1,3]
                  w4 = w[1,4]
            }
            weights = rbind(weights,c(w1,w2,w3,w4))
            temp = (w1*diff(s1) + w2*diff(s2) +w3*diff(s3)+w4*diff(s4))
            final = min(opt_to+n,2309)
            p = rbind(p,temp[(opt_to+1):final,])
            
      }
            
      return(list(eq=cumsum(p),w=weights))
}
get_sig_sma_opt<-function(para,data,threshold=0,from=1,to=500){
      mu = c(rep(0,29),para$mu)
      if(is.na(to)) to = 2309
      n = length(mu)
      thres = seq(0,1,by=0.1)
      dy = OpCl(data)
      n.grid = seq(3,60)
      mat_b = matrix(0,nrow = length(thres),ncol=length(n.grid))
      mat_s = matrix(0,nrow = length(thres),ncol=length(n.grid))
      for(j in 1:length(thres))
            for(i in 1:length(n.grid)){
                  dir_mu_b = rep(0,n)
                  dir_mu_s = rep(0,n)
                  dir_mu_b[mu>thres[j]] = 1
                  dir_mu_s[mu<(-thres[j])] = -1
                  dir_mu_b = Lag(dir_mu_b,1)
                  dir_mu_s = Lag(dir_mu_s,1)
                  dir_mu_b[1] = 0
                  dir_mu_s[1] = 0
                  sma = SMA(Op(data),n=n.grid[i])
                  sma_dir_b = rep(0,n)
                  sma_dir_s = rep(0,n)
                  sma_dir_b[sma<Op(data)] = 1
                  sma_dir_s[sma>Op(data)] = -1
                  call_b = rep(0,n)
                  call_s = rep(0,n)
                  call_b[(dir_mu_b+sma_dir_b)>1]=1
                  call_s[(dir_mu_s+sma_dir_s)<(-1)]=-1
                  trans_b = get_trans(call_b)
                  trans_s = get_trans(call_s)
                  eq = cbind(cumsum(call_b*dy-trans_b),cumsum(call_s*dy-trans_s),cumsum(call_b*dy)+cumsum(call_s*dy)-trans_b-trans_s,cumsum(dy))[from:to]
                  names(eq) = c('strat_buy','strat_sell','port','bnh')
                  stat = apply(apply(eq,2,diff),2,Return.annualized)/apply(apply(eq,2,diff),2,maxDrawdown)
                  mat_b[j,i] = stat[1]
                  mat_s[j,i] = stat[2]
      }
      return(list(buy=mat_b,sell=mat_s))
}
get_sig_sma<-function(para,data,threshold=0,from=1,to=NA,lb=9){
      mu = c(rep(0,29),para$mu)
      if(is.na(to)) to = 2309
      n = length(mu)
      dir_mu_b = rep(0,n)
      dir_mu_s = rep(0,n)
      dir_mu_b[mu>threshold] = 1
      dir_mu_s[mu<(-threshold)] = -1
      dir_mu_b = Lag(dir_mu_b,1)
      dir_mu_s = Lag(dir_mu_s,1)
      dir_mu_b[1] = 0
      dir_mu_s[1] = 0
      sma = SMA(Op(data),n=lb)
      sma_dir_b = rep(0,n)
      sma_dir_s = rep(0,n)
      sma_dir_b[sma<Op(data)] = 1
      sma_dir_s[sma>Op(data)] = -1
      call_b = rep(0,n)
      call_s = rep(0,n)
      call_b[(dir_mu_b+sma_dir_b)>1]=1
      call_s[(dir_mu_s+sma_dir_s)<(-1)]=-1
      eq = cbind(cumsum(call_b*dy),cumsum(call_s*dy),cumsum(call_b*dy)+cumsum(call_s*dy),cumsum(dy))[from:to]
      names(eq) = c('strat_buy','strat_sell','port','bnh')
      plot(as.zoo(eq),plot.type='single',col=c('blue','red','black','magenta'),main='Strat Returns vs BnH Returns',ylab='Cumulative Eq',xlab='Time')
      legend('topleft',col=c('blue','red','black','magenta'),lwd=c(2,2,2),lty=c(1,1,1),c('Buy','Sell','Comb','BnH'))
      print(apply(apply(eq,2,diff),2,Return.annualized)/apply(apply(eq,2,diff),2,maxDrawdown))
      return(eq)
}

get_sig_ema_opt<-function(para,data,threshold=0,from=1,to=500){
      mu = c(rep(0,29),para$mu)
      if(is.na(to)) to = 2309
      n = length(mu)
      thres = seq(0,1,by=0.1)
      dy = OpCl(data)
      n.grid = seq(3,30)
      mat_b = matrix(0,nrow = length(thres),ncol=length(n.grid))
      mat_s = matrix(0,nrow = length(thres),ncol=length(n.grid))
      for(j in 1:length(thres))
            for(i in 1:length(n.grid)){
                  dir_mu_b = rep(0,n)
                  dir_mu_s = rep(0,n)
                  dir_mu_b[mu>thres[j]] = 1
                  dir_mu_s[mu<(-thres[j])] = -1
                  dir_mu_b = Lag(dir_mu_b,1)
                  dir_mu_s = Lag(dir_mu_s,1)
                  dir_mu_b[1] = 0
                  dir_mu_s[1] = 0
                  sma = EMA(Op(data),n=n.grid[i])
                  sma_dir_b = rep(0,n)
                  sma_dir_s = rep(0,n)
                  sma_dir_b[sma>Op(data)] = 1
                  sma_dir_s[sma<Op(data)] = -1
                  call_b = rep(0,n)
                  call_s = rep(0,n)
                  call_b[(dir_mu_b+sma_dir_b)>1]=1
                  call_s[(dir_mu_s+sma_dir_s)<(-1)]=-1
                  eq = cbind(cumsum(call_b*dy),cumsum(call_s*dy),cumsum(call_b*dy)+cumsum(call_s*dy),cumsum(dy))[from:to]
                  names(eq) = c('strat_buy','strat_sell','port','bnh')
                  stats = apply(apply(eq,2,diff),2,Return.annualized)/apply(apply(eq,2,diff),2,maxDrawdown)
                  mat_b[j,i] = stats[1]
                  mat_s[j,i] = stats[2]
            }
      return(list(buy=mat_b,sell=mat_s))
}


get_sig_ema<-function(para,data,threshold=0,from=1,to=NA,lb=5){
      mu = c(rep(0,29),para$mu)
      if(is.na(to)) to = 2309
      n = length(mu)
      dir_mu_b = rep(0,n)
      dir_mu_s = rep(0,n)
      dir_mu_b[mu>threshold] = 1
      dir_mu_s[mu<(-threshold)] = -1
      dir_mu_b = Lag(dir_mu_b,1)
      dir_mu_s = Lag(dir_mu_s,1)
      dir_mu_b[1] = 0
      dir_mu_s[1] = 0
      dy = OpCl(data)
      sma = EMA(Op(data),n=lb)
      sma_dir_b = rep(0,n)
      sma_dir_s = rep(0,n)
      sma_dir_b[sma>Op(data)] = 1
      sma_dir_s[sma<Op(data)] = -1
      call_b = rep(0,n)
      call_s = rep(0,n)
      call_b[(dir_mu_b+sma_dir_b)>1]=1
      call_s[(dir_mu_s+sma_dir_s)<(-1)]=-1
      eq = cbind(cumsum(call_b*dy),cumsum(call_s*dy),cumsum(call_b*dy)+cumsum(call_s*dy),cumsum(dy))[from:to]
      names(eq) = c('strat_buy','strat_sell','port','bnh')
      plot(as.zoo(eq),plot.type='single',col=c('red','blue','black','magenta'),main='Strat Returns vs BnH Returns',ylab='Cumulative Eq',xlab='Time')
      legend('topleft',col=c('red','blue','black','magenta'),lwd=c(2,2,2,2),lty=c(1,1,1,1),c('Strat_buy','Strat_sell','Comb','BnH'))
      print(apply(apply(eq,2,diff),2,Return.annualized)/apply(apply(eq,2,diff),2,maxDrawdown))
      return(eq)
}
get_sig_mo_opt<-function(para,data,lb_y = 4,from=1,to=500){
      mu = c(rep(0,29),para$mu)
      if(is.na(to)) to = 2309
      n = length(mu)
      thres = seq(0,1,by=0.1)
      y_thres = seq(0,0.05,by=0.01)
      score = seq(0,lb_y)
      dy = OpCl(data)
      output =c()
      for(i in thres)
            for(j in y_thres)
                  for(w in score){
                        dir_mu_b = rep(0,n)
                        dir_mu_s = rep(0,n)
                        dir_mu_b[mu>i] = 1
                        dir_mu_s[mu<(-i)] = -1
                        dir_mu_b = Lag(dir_mu_b,1)
                        dir_mu_s = Lag(dir_mu_s,1)
                        dir_mu_b[1] = 0
                        dir_mu_s[1] = 0
                        dir = rep(0,n)
                        dir[dy>j] = 1
                        dir[dy<(-j)] = -1
                        dir = Lag(dir,1)
                        dir[1] = 0
                        dir = filter(dir,rep(1,lb_y), sides=1)
                        call_b = rep(0,n)
                        call_s = rep(0,n)
                        call_b[(dir_mu_b+dir)>w] = 1
                        call_s[(dir_mu_s+dir)<(-w)] = -1
                        buy = call_b*dy
                        sell = call_s*dy
                        comb = cbind(buy,sell,buy+sell)
                        temp = apply(apply(comb,2,diff),2,Return.annualized)/apply(apply(comb,2,diff),2,maxDrawdown)
                        output = rbind(output,c(i,j,w,temp))
                              
                        
                  }
      
      
      return(output)
}
get_sig_mo<-function(para,data,threshold=0,y_thres=0,lb_y=3,score=2,from=1,to=NA){
      mu = c(rep(0,29),para$mu)
      if(is.na(to)) to = 2309
      n = length(mu)
      dir_mu_b = rep(0,n)
      dir_mu_s = rep(0,n)
      dir_mu_b[mu>threshold] = 1
      dir_mu_s[mu<(-threshold)] = -1
      dir_mu_b = Lag(dir_mu_b,1)
      dir_mu_s = Lag(dir_mu_s,1)
      dir_mu_b[1] = 0
      dir_mu_s[1] = 0
      dy = OpCl(data)
      dir = rep(0,n)
      dir[dy>y_thres] = 1
      dir[dy<(-y_thres)] = -1
      dir = Lag(dir,1)
      dir[1] = 0
      dir = filter(dir,rep(1,lb_y), sides=1)
      call_b = rep(0,n)
      call_s = rep(0,n)
      call_b[(dir_mu_b+dir)>score]=1
      call_s[(dir_mu_s+dir)<(-score)]=-1
      eq = cbind(cumsum(call_b*dy),cumsum(call_s*dy),cumsum(call_b*dy+call_s*dy),cumsum(dy))[from:to]
      names(eq) = c('strat_buy','strat_sell','port','bnh')
      plot(as.zoo(eq),plot.type='single',col=c('red','blue','black','magenta'),main='Strat Returns vs BnH Returns',ylab='Cumulative Eq',xlab='Time')
      legend('topleft',col=c('red','blue','black','magenta'),lwd=c(2,2,2,2),lty=c(1,1,1,1),c('Strat_buy','Strat_sell','Comb','BnH'))
      print(apply(apply(eq,2,diff),2,Return.annualized)/apply(apply(eq,2,diff),2,maxDrawdown))
      return(eq)
}

get_sig_naive<-function(para,data,threshold=0,from=1,to=NA){
      mu = c(rep(0,29),para$mu)
      if(is.na(to)) to = 2309
      n = length(mu)
      dir_mu_b = rep(0,n)
      dir_mu_s = rep(0,n)
      dir_mu_b[mu>threshold] = 1
      dir_mu_s[mu<(-threshold)] = -1
      dy = OpCl(data)
      dir_mu_b = Lag(dir_mu_b,1)
      dir_mu_s = Lag(dir_mu_s,1)
      dir_mu_b[1] = 0
      dir_mu_s[1] = 0
      trans_b = get_trans(call_b)
      trans_s = get_trans(call_s)
      eq = cbind(cumsum(call_b*dy-trans_b),cumsum(call_s*dy-trans_s),cumsum(call_b*dy)+cumsum(call_s*dy)-trans_b-trans_s,cumsum(dy))[from:to]
      names(eq) = c('strat_buy','strat_sell','port','bnh')
      plot(as.zoo(eq),plot.type='single',col=c('red','blue','black','magenta'),main='Strat Returns vs BnH Returns',ylab='Cumulative Eq',xlab='Time')
      legend('topleft',col=c('red','blue','black','magenta'),lwd=c(2,2,2,2),lty=c(1,1,1,1),c('Strat_buy','Strat_sell','Comb','BnH'))
      print(apply(apply(eq,2,diff),2,Return.annualized)/apply(apply(eq,2,diff),2,maxDrawdown))
      return(eq)
}


