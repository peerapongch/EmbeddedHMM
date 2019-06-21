load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')
source('../../samplers/model1_mhmcmc.R')

### MHMCMC with autoregressive update ###
N <- 2800000
es <- c(0.2,0.8)
seed <- 1

system.time(
  mcmc_out <- mcmcGaussianSSM(N,es,ssm_poisson,obs='Poisson',seed=seed) 
)

save(mcmc_out,file='./data/test1_poisson1_mcmc_N2800000_seed1.RData')
# user   system  elapsed 
# 14803.47    10.53 15058.29 
# tps = 0.0537

# adjustment for time 
mcmc_out$X_sample <- mcmc_out$X_sample[1:150000,,]

### compute posterior mean ###
mci_mu <- matrix(0,nrow=T,ncol=dim)
for(t in 1:T){
  mci_mu[t,] <- mean(mcmc_out$X_sample[,t,])
}

### plot samples ###
source('./function_plot_mcmc.R')
par(mfrow=c(1,1))
for(i in 1:dim){
  plot_mcmc_mu(i,ssm_poisson,mcmc_out,mci_mu,plot.samples=TRUE,interval=1000) 
  # plot(ssm_poisson$Y[,i],type='l')
}

### visualise autocorrelation ###
acf(mcmc_out$X_sample[,100,1],lag.max=100)
# use lag.max=50
source('../../evaluation/autocorrelation_time.R')
actime_out <- ACTime(list(mcmc_out),T,dim,60,0.0537)
plot(actime_out$ac_time[1,,1],type='l',main='dim 1')
plot(actime_out$ac_time[1,,2],type='l',main='dim 2')
plot(actime_out$ac_time[1,,3],type='l',main='dim 3')
