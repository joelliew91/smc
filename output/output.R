
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
