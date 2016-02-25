source('~/Desktop/FYP/sim_svvg.r')


for(i in 1:100){
      data = sim_svvg(150)
      write.csv(data$y,sprintf('~/Desktop/smc/data/data%d%s',i,'.csv'))
}