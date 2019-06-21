load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')
source('../../samplers/model1_mhmcmc.R')

es <- c(0.2,0.8)
N <- 10
seed <- 1
system.time(
  mcmc_out <- mcmcGaussianSSM(N,es,ssm_poisson,seed=seed)
)
# user  system elapsed 
# 46.25    0.10   47.36 
# tpi = 0.0047
tps <- 47.36/N*10 #

### compute posterior mean ###
mci_mu <- matrix(0,nrow=T,ncol=dim)
for(t in 1:T){
  mci_mu[t,] <- mean(mcmc_out$X_sample[,t,])
}

### plot samples ###
source('./function_plot_mcmc.R')
par(mfrow=c(1,1))
for(i in 1:dim){
  plot_mcmc_mu(i,ssm_poisson,mcmc_out,mci_mu,plot.samples=TRUE,interval=200) 
  # plot(ssm_poisson$Y[,i],type='l')
}

### visualise autocorrelation ###
acf(mcmc_out$X_sample[,3,2],lag.max=20)

### traceplot ### 
library(coda)
X_mcmc <- as.mcmc(mcmc_out$X_sample[,3,])
plot(X_mcmc)

# use lag.max=20
source('../../evaluation/autocorrelation_time.R')
actime_out <- ACTime(list(mcmc_out),T,dim,100,1)
plot(actime_out$ac_time[1,,1],type='l',main='dim 1')
plot(actime_out$ac_time[1,,2],type='l',main='dim 2')
plot(actime_out$ac_time[1,,3],type='l',main='dim 3')
