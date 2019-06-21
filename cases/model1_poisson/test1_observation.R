load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')

for(i in 1:dim){
  plot(ssm_poisson$Y[,i],type='l',
       main=paste('Dim',i,sep=''),ylab='')
  lines(exp(c[i]+delta[i]*ssm_poisson$X[,i]),col='red')
  legend(0,10,legend=c('Observation','Prior mean'),
         col=c('black','red'), lty=c(1,1),cex=0.8)
}
