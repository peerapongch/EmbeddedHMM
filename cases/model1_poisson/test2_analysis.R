load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')

load(file='./data/test2_poisson1_pgmet_L20_N17000_seed1.RData')
load(file='./data/test2_poisson1_pgbs_L20_N20000_seed1.RData')

### adjustment ### 
pgmet_out$ X_sample <- pgmet_out$X_sample[1:8500,,]

### calculate means ###
mcmcs <- list(pgbs_out,pgmet_out)
mci_mus <- array(0,dim=c(3,T,dim))
for(i in 1:length(mcmcs)){
  for(t in 1:T){
    mci_mus[i,t,] <- mean(mcmcs[[i]]$X_sample[,t,])
  }
}

### old mus ### 
rm(pgbs_out); rm(pgmet_out)
load('./data/test1_poisson1_pgbs_L250_N20000_seed1.RData')
load('./data/test1_poisson1_pgmet_L250_N20000_seed1.RData')
pgmet_out$X_sample <- pgmet_out$X_sample[1:64000,,]
mcmcs_old <- list(pgbs_out,pgmet_out)
mci_mus_old <- array(0,dim=c(3,T,dim))
for(i in 1:length(mcmcs_old)){
  for(t in 1:T){
    mci_mus_old[i,t,] <- mean(mcmcs_old[[i]]$X_sample[,t,])
  }
}


### compare mus ###
for(i in 1:dim){
  plot(ssm_poisson$X[,i],type='l',ylim=c(-6,10),
       main ='Comparison of approximated smoothing mean',
       ylab = paste('x',i,sep='_'))
  lines(mci_mus[1,,i],col='blue')
  lines(mci_mus[2,,i],col='green')
  lines(mci_mus_old[1,,i],col='red')
  lines(mci_mus_old[2,,i],col='orange')
  legend(0,10,legend=c('True latent state','PGBS (L=20)','PGBS+Metropolis (L=20)','PGBS (L=250)','PGBS+Metropolis (L=250)'),
         col=c('black','blue','green','red','orange'), lty=c(1,1),cex=0.8)
}

### pgbs acf ### 
acf(pgbs_out$X_sample[,3,1],lag.max=10)
acf(pgbs_out$X_sample[,3,2],lag.max=10)
acf(pgbs_out$X_sample[,3,3],lag.max=10)

### pgmet acf ### 
acf(pgmet_out$X_sample[,3,1],lag.max=20)
acf(pgmet_out$X_sample[,3,2],lag.max=20)
acf(pgmet_out$X_sample[,3,3],lag.max=20)

### unadjusted actime ###
source('../../evaluation/autocorrelation_time.R')
actimes <- list()
tps <- c(1,1)
actimes[[1]] <- ACTime(list(mcmcs[[1]]),T,dim,20,tps[1])
actimes[[2]] <- ACTime(list(mcmcs[[2]]),T,dim,20,tps[1])
actimes[[3]] <- ACTime(list(mcmcs_old[[1]]),T,dim,20,tps[1])
actimes[[4]] <- ACTime(list(mcmcs_old[[2]]),T,dim,20,tps[1])
for(i in 1:dim){
  plot(actimes[[1]]$ac_time[1,,i],type='l',ylim=c(0,10),
       main='Unadjusted autocorrelation times',ylab=paste('x_',i,sep=''),col='blue')
  lines(actimes[[2]]$ac_time[1,,i],col='green')
  # plot(actimes[[2]]$ac_time[1,,i],type='l',ylim=c(0,4),
  # main='Unadjusted autocorrelation times',ylab=paste('x_',i,sep=''),col='red')
  lines(actimes[[3]]$ac_time[1,,i],col='red')
  lines(actimes[[4]]$ac_time[1,,i],col='orange')
  legend(0,10,legend=c('PGBS (L=20)','PGBS+Metropolis (L=20)','PGBS (L=250)','PGBS+Metropolis (L=250)'),
         col=c('blue','green','red','orange'), lty=c(1,1),cex=0.8)
}

### Adjusted actime ###
source('../../evaluation/autocorrelation_time.R')
actimes <- list()
tps <- c(0.05,0.19,0.12)
tps <- c(1,1,1)
actimes[[1]] <- ACTime(list(mcmc_out),T,dim,60,tps[1])
actimes[[2]] <- ACTime(list(pgbs_out),T,dim,10,tps[2])
actimes[[3]] <- ACTime(list(pgmet_out),T,dim,10,tps[3])
for(i in 1:dim){
  plot(actimes[[1]]$ac_time[1,,i],type='l',ylim=c(0,1),
       main='Unadjusted autocorrelation times',ylab=paste('x_',i,sep=''))
  lines(actimes[[2]]$ac_time[1,,i],col='red')
  # plot(actimes[[2]]$ac_time[1,,i],type='l',ylim=c(0,4),
  # main='Unadjusted autocorrelation times',ylab=paste('x_',i,sep=''),col='red')
  lines(actimes[[3]]$ac_time[1,,i],col='blue')
  # legend(0,4,legend=c('PGBS','PGBS+Metropolis'),
  # col=c('red','blue'), lty=c(1,1),cex=0.8)
  legend(150,1,legend=c('Metropolis','PGBS','PGBS+Metropolis'),
         col=c('black','red','blue'), lty=c(1,1),cex=0.8)
}

