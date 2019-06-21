source('../../samplers/model1_pgmet.R')
load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')

L <- 20
N <- 17000
seed <- 1
es <- c(0.2,0.8)

system.time(
  pgmet_out <- pgmetModel1(ssm_poisson,N,L,es,seed=seed)
)
# user  system elapsed 
# 1997.70    2.82 2029.44
tps <- 0.0298

save(pgmet_out,file='./data/test2_poisson1_pgmet_L20_N17000_seed1.RData')

### adjustment ### 
pgmet_out$ X_sample <- pgmet_out$X_sample[1:8500,,]

### compute posterior mean ###
mci_mu <- matrix(0,nrow=T,ncol=dim)
for(t in 1:T){
  mci_mu[t,] <- mean(pgmet_out$X_sample[,t,])
}

### plot samples ###
source('./function_plot_mcmc.R')
par(mfrow=c(1,1))
for(i in 1:dim){
  plot_mcmc_mu(i,ssm_poisson,pgmet_out,mci_mu,plot.samples=TRUE,interval=200) 
  # plot(ssm_poisson$Y[,i],type='l')
}

### visualise autocorrelation ###
acf(pgmet_out$X_sample[,3,3],lag.max=20)

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

