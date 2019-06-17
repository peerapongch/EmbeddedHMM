# assuming the current working directory is at current file path 
source('../../models/model1_generate.R')
source('../../samplers/model1_mhmcmc.R')
source('../../samplers/model1_kalman.R')

# model 1 specification and simulation
T <- 250
dim <- 3
mu_init <- rep(0,dim)
rho <- 0.7; phi <- 0.9; phis <- rep(phi,dim); v <- 1/sqrt((1-phis^2)); sigma_init <- v %*% t(v)
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
delta <- 0.6; dim <- dim; deltas <- rep(delta,dim)

F <- diag(phis); G <- t(chol(sigma)); R <- diag(1,dim); Q <- diag(1,dim); H <- diag(deltas)
ssm <- generateGaussianGaussianSSM(T,dim,mu_init,sigma_init,F,G,Q,H,R)
save(ssm,file='./data/ssm_model1_gaussian_1.RData')

load('./data/ssm_model1_gaussian_1.RData')

### Kalman filtering ### 
kfks <- KFKS(ssm)

### MHMCMC with autoregressive proposal ### 
N <- 1000; e <- c(0.2,0.8)
system.time(
  mcmc <- mcmcGaussianSSM(N,es,ssm)  
)
X_mcmc <- mcmc$X_mcmc

### plot diagnostics ### 
mci_mu <- matrix(0,nrow=ssm$T,ncol=ssm$dim)
for(j in 1:ssm$dim){
  for(i in 1:ssm$T){
    mci_mu[i,j] <- mean(X_mcmc[,i,j])
  }
}
plot_compare_gaussian <- function(d,ssm,kfks,mcmc,mci_mu){
  X_mcmc <- mcmc$X_mcmc
  plot(X_mcmc[1,,d],type='l',col='grey',ylim=c(-6,7),ylab=paste('X_',d,sep=''),xlab='t',main='Comparison of Posterior smoothing means')
  for(i in 2:N){
    lines(X_mcmc[i,,d],col='grey')
  }
  lines(mci_mu[,d],col='black')
  lines(ssm$X[,d],type='l',col='blue')
  lines(kfks$mu_smoothing[,d],col='red')
  legend(0,7,legend=c('MHMCMC post. mean','Exact post. mean','Latent state (truth)'),
         col=c('black','red','blue'), lty=c(1,1,1),cex=0.8)
}
### dim 
for(i in 1:dim){
  plot_compare_gaussian(i,ssm,kfks,mcmc,mci_mu)
}
