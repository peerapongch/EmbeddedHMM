load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')
source('../../samplers/model1_mhmcmc.R')

### MHMCMC with autoregressive update ###
N <- 10000; es <- c(0.2,0.8)
system.time(
  mhmcmc_out <- mhmcmcModel1(ssm_poisson,N,es) 
)

#  save(mcmc,file='./data/model1_poisson1_mcmc_N10e6.RData')

X_mcmc <- mcmc$X_sample
### MCI mu ###
mci_mu <- matrix(0,nrow=ssm_poisson$T,ncol=ssm_poisson$dim)
for(j in 1:ssm_poisson$dim){
  for(i in 1:ssm_poisson$T){
    mci_mu[i,j] <- mean(X_mcmc[,i,j])
  }
}

### diagnostic (compare with actual latent states) ###
source('./function_plot_mcmc.R')

par(mfrow=c(1,1))
for(i in 1:dim){
  plot_mcmc_mu(i,ssm_poisson,mcmc,mci_mu,plot.samples=TRUE,interval=100) 
  # plot(ssm_poisson$Y[,i],type='l')
}

# correlation? visualise the samples 
library(MCMCpack)
x <- X_mcmc[,1:3,1]
length(unique(x[,1]))/length(x[,1])
length(unique(x[,2]))/length(x[,2])
length(unique(x[,3]))/length(x[,3])
x <- as.mcmc(x)
plot(x)
