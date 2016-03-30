sim_svvg<-function(years,delta=1){
    y = c(1)                                        ## Y0 set at 0
    vol = c(1)                                      ## Volatility at t=0 is set at 1
    mu = 0.05
    k = 0.015
    theta = 0.8
    sigma.v = 0.1
    rho = -0.4
    gam = -0.01
    sigma.y = 0.4
    v = 3.0
    sigma = 0
    E = matrix(c(1,rho*sigma.v,rho*sigma.v,sigma.v^2),2,2)
    epsilon = mvrnorm(years*250,mu=c(0,0),Sigma=E)  ## Generate 5000 moves ie 20yrs in the
    ## BM of log price and volatility
    epsilon.j = rnorm(years*250)                    ## Generate the Jump outcomes
    Gamma_path = rgamma(years*250,delta/v,v)
    J = c()
    
    for(i in Gamma_path){
        temp = rnorm(1,mean=gam*i,sd=sigma.v*sqrt(i))
        J = c(J,temp)
    }
    
    for(i in 1:(years*250)){
        Y = y[i] + mu*delta + sqrt(vol[i]*delta)*epsilon[i,1] + J[i]
        V = vol[i] + k*(theta-vol[i])*delta + sigma.v*sqrt(vol[i]*delta)*epsilon[i,2]
        y = rbind(y,Y)
        vol = rbind(vol,V)
    }
    par(mfrow=c(2,1))
    plot(y,type='l',ylab='Log Stock Price',xlab = 'Time',main='SVVG')
    plot(exp(y),type='l',ylab='Stock Price',xlab = 'Time',main = 'Stock Price using SVVG')
    return(list(y=y,vol=vol))
}






