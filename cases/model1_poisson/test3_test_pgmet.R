load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')
source('../../samplers/model1_pgmet.R')

L <- 250
N <- 100
seed <- 1

system.time(
  pgmet_out <- pgmetModel1(ssm_poisson,N,L,seed)
)
# user  system elapsed 
# 108.37    9.23   64.89  
# tpi = 0.65
tps <- 64.89/N/2 # 0.30

# save(pgmet_out,file='./data/test2_poisson1_pgbs_L20_N20000_seed1.RData')

### compute posterior mean ###
mci_mu <- matrix(0,nrow=T,ncol=dim)
for(t in 1:T){
  mci_mu[t,] <- mean(pgmet_out$X_sample[,t,])
}

### plot samples ###
source('./function_plot_mcmc.R')
par(mfrow=c(1,1))
for(i in 1:dim){
  plot_mcmc_mu(i,ssm_poisson,pgmet_out,mci_mu,plot.samples=TRUE,interval=10) 
  # plot(ssm_poisson$Y[,i],type='l')
}

### visualise autocorrelation ###
acf(pgmet_out$X_sample[,3,2],lag.max=200)

### traceplot ### 
library(coda)
X_mcmc <- as.mcmc(pgmet_out$X_sample[,3,])
plot(X_mcmc)

# use lag.max=20
source('../../evaluation/autocorrelation_time.R')
actime_out <- ACTime(list(pgmet_out),T,dim,20,1)
plot(actime_out$ac_time[1,,1],type='l',main='dim 1')
plot(actime_out$ac_time[1,,2],type='l',main='dim 2')
plot(actime_out$ac_time[1,,3],type='l',main='dim 3')
