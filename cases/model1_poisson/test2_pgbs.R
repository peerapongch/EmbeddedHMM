load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')
source('../../samplers/model1_pgbs.R')

L <- 20
N <- 20000
seed <- 1

system.time(
  pgbs_out <- pgbsModel1(ssm_poisson,N,L,seed)
)
# user  system elapsed 
# 1043.19    0.53 1049.18 
tps <- 0.0263

save(pgbs_out,file='./data/test2_poisson1_pgbs_L20_N20000_seed1.RData')

### compute posterior mean ###
mci_mu <- matrix(0,nrow=T,ncol=dim)
for(t in 1:T){
  mci_mu[t,] <- mean(pgbs_out$X_sample[,t,])
}

### plot samples ###
source('./function_plot_mcmc.R')
par(mfrow=c(1,1))
for(i in 1:dim){
  plot_mcmc_mu(i,ssm_poisson,pgbs_out,mci_mu,plot.samples=TRUE,interval=200) 
  # plot(ssm_poisson$Y[,i],type='l')
}

### visualise autocorrelation ###
acf(pgbs_out$X_sample[,3,2],lag.max=20)

### traceplot ### 
library(coda)
X_mcmc <- as.mcmc(pgbs_out$X_sample[,3,])
plot(X_mcmc)

# use lag.max=20
source('../../evaluation/autocorrelation_time.R')
actime_out <- ACTime(list(pgbs_out),T,dim,20,1)
plot(actime_out$ac_time[1,,1],type='l',main='dim 1')
plot(actime_out$ac_time[1,,2],type='l',main='dim 2')
plot(actime_out$ac_time[1,,3],type='l',main='dim 3')

