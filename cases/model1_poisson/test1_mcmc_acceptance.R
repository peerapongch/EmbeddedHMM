load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')
source('../../samplers/model1_mhmcmc.R')

### MHMCMC with autoregressive update ###
N <- 10000
es <- c(0.2,0.8)
seed <- 1

system.time(
  mcmc_out <- mcmcGaussianSSM(N,es,ssm_poisson,obs='Poisson',seed=seed,thin.factor=1) 
)

acceptance_rate <- matrix(0,nrow=T,ncol=dim)
for(t in 1:T){
  for(j in  1:dim){
    acceptance_rate[t,j] <- length(unique(mcmc_out$X_sample[,t,j]))/length(mcmc_out$X_sample[,t,j])
  }
}
plot(acceptance_rate[,1],type='l',main='acceptance rate')


library(coda)
arate <- rep(0,T)
for(t in 1:T){
  mcmc <- as.mcmc(mcmc_out$X_sample[,t,1])
  arate[t] <- 1-rejectionRate(mcmc)[[1]]
}
plot(arate,type='l')
