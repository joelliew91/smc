
get_para_all<-function(){
      o = read.csv(sprintf('output/output%d%s',1,'.csv'),header=T)
      for(i in 2:100){
            o = rbind(o,read.csv(sprintf('output/output%d%s',i,'.csv'),header=T))
      }
      return(o)
}

get_para<-function(){
    para = as.data.frame(matrix(0,nrow=100,ncol=7))
    names(para) = c('mu','gam','sigma_j','sigma_v','rho','k','v')
    for(i in 1:100){
        o = read.csv(sprintf('~/Desktop/smc/output/output%d%s',i,'.csv'),header=T)
        for(j in 1:7){
            para[i,j] = sum(o[,j]*o$norm_w)
        }
    }
    return(para)
}

get_hist<-function(para){
      p = c(0.085,0.0,0.4,0.14,-0.1,1.6,0.155)
      par(mfrow=c(3,3))
      breaks = 100
      if(length(para[,1])<200)
            breaks = 30
      for(i in 1:7){

#            if(i==4){
#                  hist(para[,i],breaks=breaks,main = sprintf("Histogram of %s",names(para)[i]),xlim=c(0,0.5))  
#            }
#            else{
#                  hist(para[,i],breaks=breaks,main = sprintf("Histogram of %s",names(para)[i]))  
#            }     
            plot(density(para[,i]),main = sprintf("Histogram of %s",names(para)[i]))



            abline(v=quantile(para[,i],0.975),col='blue',lwd=3,lty=3)
            abline(v=quantile(para[,i],0.025),col='blue',lwd=3,lty=3)
            abline(v=p[i],col='red',lwd=3)
      }

      
}