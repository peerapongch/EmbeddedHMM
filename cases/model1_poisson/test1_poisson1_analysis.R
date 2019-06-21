load('./data/test1_poisson1_mcmc_N2800000_seed1.RData')
load('./data/test1_poisson1_pgbs_L250_N20000_seed1.RData')
load('./data/test1_poisson1_pgmet_L250_N20000_seed1.RData')
load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')

# check 
dim(pgbs_out$X_sample)
dim(pgmet_out$X_sample)
dim(mcmc_out$X_sample)

### adjust ###
# pgbs_out$X_sample <- pgbs_out$X_sample[1:10,,]
mcmc_out$X_sample <- mcmc_out$X_sample[1:150000,,]
pgmet_out$X_sample <- pgmet_out$X_sample[1:64000,,]

### calculate means ###
mcmcs <- list(mcmc_out,pgbs_out,pgmet_out)
mci_mus <- array(0,dim=c(3,T,dim))
for(i in 1:3){
  for(t in 1:T){
    mci_mus[i,t,] <- mean(mcmcs[[i]]$X_sample[,t,])
  }
}

mci_mus[1,1,]
mci_mus[2,1,]
mci_mus[3,1,]

### compare mus ###
for(i in 1:dim){
  plot(ssm_poisson$X[,i],type='l',ylim=c(-6,10),
       main ='Comparison of approximated smoothing mean',
       ylab = paste('x',i,sep='_'))
  lines(mci_mus[1,,i],col='red')
  lines(mci_mus[2,,i],col='blue')
  lines(mci_mus[3,,i],col='green')
  legend(0,10,legend=c('True latent state','Metropolis','PGBS','PGBS+Metropolis'),
         col=c('black','red','blue','green'), lty=c(1,1),cex=0.8)
}
# literal chills 

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
tps <- c(1,1,1)
actimes[[1]] <- ACTime(list(mcmc_out),T,dim,60,tps[1])
actimes[[2]] <- ACTime(list(pgbs_out),T,dim,10,tps[2])
actimes[[3]] <- ACTime(list(pgmet_out),T,dim,10,tps[3])
for(i in 1:dim){
  plot(actimes[[1]]$ac_time[1,,i],type='l',ylim=c(0,20),
       main='Unadjusted autocorrelation times',ylab=paste('x_',i,sep=''))
  lines(actimes[[2]]$ac_time[1,,i],col='red')
  # plot(actimes[[2]]$ac_time[1,,i],type='l',ylim=c(0,4),
       # main='Unadjusted autocorrelation times',ylab=paste('x_',i,sep=''),col='red')
  lines(actimes[[3]]$ac_time[1,,i],col='blue')
  # legend(0,4,legend=c('PGBS','PGBS+Metropolis'),
         # col=c('red','blue'), lty=c(1,1),cex=0.8)
  legend(0,20,legend=c('Metropolis','PGBS','PGBS+Metropolis'),
         col=c('black','red','blue'), lty=c(1,1),cex=0.8)
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

