sim_svvg<-function(years,delta=1/250){
      y = c(1)                                        ## Y0 set at 0
      vol = c(0.1)                                      ## Volatility at t=0 is set at 1
      mu = 0.085
      k = 1.6
      v_p = 0.1
      sigma.v = 0.14
      rho = -0.1
      gam = 0.0
      sigma.j = 0.4
      
      Gamma_path = rgamma(years,1,1)
      Z = c()
      
      for(i in Gamma_path){
            temp = -gam*i+sigma.j*sqrt(i)*rnorm(1,mean = 0,sd=1)
            Z = c(Z,temp)
      }
      
      for(i in 1:(years)){
            u1 = rnorm(1)
            V = abs(vol[i] + k*(v_p-vol[i])*delta + sigma.v*sqrt(vol[i]*delta)*u1)
            Y = y[i] + mu*delta + sqrt(vol[i]*delta)*(rho*u1+sqrt(1-rho*rho)*rnorm(1)) + Z[i]
            
            y = rbind(y,Y)
            vol = rbind(vol,V)
      }
      par(mfrow=c(3,1))
      plot(y,type='l',ylab='Log Stock Price',xlab = 'Time',main='SVVG')
      plot(vol,type='l',ylab='Vol',xlab = 'Time',main = 'Vol using SVVG')
      plot(Z,type='l',ylab='Variance Gamma',xlab='Time',main='VG')
      return(list(y=y,vol=vol))
}






