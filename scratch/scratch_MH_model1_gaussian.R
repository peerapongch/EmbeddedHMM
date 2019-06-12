install.packages('mvtnorm')
library(mvtnorm)

pre_mu_1 <- solve(F %*% F + solve(sigma_init)%*%sigma) %*% F
sigma_1 <- solve(F %*% (solve(sigma) %*% F) + solve(sigma_init))

pre_mu_j <- solve(F %*% F + diag(1,dim(F)[1])) %*% F
sigma_j <- solve(F %*% (solve(sigma) %*% F) + solve(sigma))

pre_mu_T <- F
sigma_T <- sigma

##### MHMCMC autoregressive ######
N <- 1000; e <- 0.5 # ranges from -1 to 1
X_mcmc <- array(0,dim=c(N,T,dim))
X_mcmc[1,,] <- mvrnorm(T,mu_init,sigma_init)
pb <- txtProgressBar(min=0,max=N,title="MH-MCMC",style=3)
for(i in 2:N){
  # print(paste('progress: ',i/N*100,'%',sep=''))
  setTxtProgressBar(pb, i)
  for(j in 1:T){
    # find the mu
    if(j==1) {
      mu_j <- pre_mu_1 %*% X_mcmc[i-1,2,]
      L <- t(chol(sigma_1))
    } else if(j==T) {
      mu_j <- pre_mu_T %*% X_mcmc[i-1,T-1,]
      L <- t(chol(sigma_T))
    } else {
      mu_j <- pre_mu_j %*% (X_mcmc[i,j-1,]+X_mcmc[i-1,j+1,])
      L <- t(chol(sigma_j))
    }

    # autoregressive update
    z <- mvrnorm(1,rep(0,dim),diag(1,dim))
    x_j <- mu_j+sqrt(1-e^2)*(X_mcmc[i-1,j,]-mu_j)+e*L%*%z

    # transition probability
    hastings <- dmvnorm(Y[j,],H %*% x_j, R)/dmvnorm(Y[j,],H %*% X_mcmc[i-1,j,], R)
    alpha <- min(1,hastings)
    U <- runif(1,0,1)
    if(alpha>U){
      X_mcmc[i,j,] <- x_j
    } else {
      X_mcmc[i,j,] <- X_mcmc[i-1,j,]
    }
  }
  close(pb)
}

mcmc_t1 <- as.mcmc(X_mcmc[,1,])

# plot(X_mcmc[1,,1],type='l')
# plot(X_mcmc[,1,1],type='l')

### plot multiple ###
# dim 1

# mcmc <- mcmcGaussianGaussianSSM(1000,0.6,ssm)

mci_mu <- rep(0,T)
for(i in 1:T){
  mci_mu[i] <- mean(X_mcmc[,i,1])
}

plot(X_mcmc[1,,1],type='l',col='grey',ylim=c(-6,7),ylab='X_1',xlab='t',main='Comparison of Posterior means')
for(i in 2:N){
  lines(X_mcmc[i,,1],col='grey')
}
lines(mci_mu,col='black')
# compare 
lines(X[,1],type='l',col='blue')
lines(mu_smoothing[,1],col='red')
legend(180,7,legend=c('MHMCMC post. mean','Exact post mean','Latent state (truth)'),
       col=c('black','red','blue'), lty=c(1,1,1),cex=0.8)

# dim 2
mci_mu <- rep(0,T)
for(i in 1:T){
  mci_mu[i] <- mean(X_mcmc[,i,2])
}

plot(X_mcmc[1,,2],type='l',col='grey',ylim=c(-8,7),ylab='X_2',xlab='t',main='Comparison of Posterior means')
for(i in 2:N){
  lines(X_mcmc[i,,2],col='grey')
}
lines(mci_mu,col='black')
# compare 
lines(X[,2],type='l',col='blue')
lines(mu_smoothing[,2],col='red')
legend(140,-4,legend=c('MHMCMC post. mean','Exact post mean','Latent state (truth)'),
       col=c('black','red','blue'), lty=c(1,1,1),cex=0.8)

# dim 3
mci_mu <- rep(0,T)
for(i in 1:T){
  mci_mu[i] <- mean(X_mcmc[,i,3])
}

plot(X_mcmc[1,,3],type='l',col='grey',ylim=c(-7,9),ylab='X_3',xlab='t',main='Comparison of Posterior means')
for(i in 2:N){
  lines(X_mcmc[i,,3],col='grey')
}
lines(mci_mu,col='black')
# compare 
lines(X[,3],type='l',col='blue')
lines(mu_smoothing[,3],col='red')
legend(50,8,legend=c('MHMCMC post. mean','Exact post mean','Latent state (truth)'),
       col=c('black','red','blue'), lty=c(1,1,1),cex=0.8)
