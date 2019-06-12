source('../../models/model1_generate.R')
source('../../samplers/model1_mhmcmc.R')

# model 1 specification and simulation: poisson
T <- 250; dim <- 3
mu_init <- rep(0,dim); rho <- 0.7; phi <- 0.9; phis <- rep(phi,dim); v <- 1/sqrt((1-phis^2)); sigma_init <- v %*% t(v)
for(i in 1:dim){
  for(j in 1:dim){
    if(i!=j){
      sigma_init[i,j] <- sigma_init[i,j] * rho
    }
  }
}
sigma <- matrix(rho,ncol=dim,nrow=dim)
for(i in 1:dim){
  sigma[i,i] <- 1
}
delta <- rep(0.6,dim); c <- rep(-0.4,dim); F <- diag(phis); G <- t(chol(sigma)); R <- diag(1,dim); Q <- diag(1,dim)
ssm_poisson <- generatePoissonGaussianSSM(T,dim,mu_init,sigma_init,F,G,Q,c,delta)
# plot(ssm_poisson$Y[,1],type='l')
save(ssm_poisson,file='./data/ssm_model1_poisson_1.RData')
# load('./data/ssm_model1_poisson_1.R')

### MHMCMC with autoregressive update ###
N <- 1000; e <- 0.5
mcmc <- mcmcPoissonGaussianSSM(N,e,ssm_poisson)
X_mcmc <- mcmc$X_mcmc

### MCI mu ###
mci_mu <- matrix(0,nrow=ssm_poisson$T,ncol=ssm_poisson$dim)
for(j in 1:ssm_poisson$dim){
  for(i in 1:ssm_poisson$T){
    mci_mu[i,j] <- mean(X_mcmc[,i,j])
  }
}

### diagnostic (compare with actual latent states) ###
plot_compare_poisson <- function(d,ssm,mcmc,mci_mu,plot.samples=FALSE,interval=1000){
  N <- mcmc$N
  thin <- round(N/interval)
  X_mcmc <- mcmc$X_mcmc
  plot(mci_mu[,d],col='black',type='l',ylim=c(-10,10),ylab=paste('X_',d,sep=''),xlab='t',
       main='Comparison of Posterior smoothing means')
  if(plot.samples){
    for(i in seq(1,N,thin)){
      lines(X_mcmc[i,,d],col='grey')
    }
    lines(mci_mu[,d],col='black')
  }
  lines(ssm$X[,d],type='l',col='blue')
  legend(0,10,legend=c('MHMCMC post. mean','Latent state (truth)'),
         col=c('black','blue'), lty=c(1,1),cex=0.8)
}

for(i in 1:dim){
  plot_compare_poisson(i,ssm_poisson,mcmc,mci_mu,plot.samples=TRUE) 
}
