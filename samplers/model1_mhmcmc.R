mcmcGaussianGaussianSSM<- function(N,e,ssm){
  # autoregressive proposal
  require(mvtnorm)
  # extract ssm values
  F <- ssm$F; T <- ssm$T; H <- ssm$H; R <- ssm$R; Q <- ssm$Q
  mu_init <- ssm$mu_init; sigma_init <- ssm$sigma_init
  sigma <- ssm$G %*% ssm$Q %*% t(ssm$G)
  dim <- ssm$dim; Y <- ssm$Y
  
  pre_mu_1 <- solve(F %*% F + solve(sigma_init)%*%sigma) %*% F
  sigma_1 <- solve(F %*% (solve(sigma) %*% F) + solve(sigma_init))
  pre_mu_j <- solve(F %*% F + diag(1,dim(F)[1])) %*% F
  sigma_j <- solve(F %*% (solve(sigma) %*% F) + solve(sigma))
  pre_mu_T <- F; sigma_T <- sigma
  
  ##### MHMCMC autoregressive ######
  # N <- 1000; e <- 0.5 # ranges from -1 to 1 
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
  }
  close(pb)
  return(list(X_mcmc=X_mcmc,N=N,e=e))
}

mcmcPoissonGaussianSSM<- function(N,e,ssm){
  # autoregressive proposal
  require(mvtnorm)
  # extract ssm values
  F <- ssm$F; T <- ssm$T; Q <- ssm$Q
  mu_init <- ssm$mu_init; sigma_init <- ssm$sigma_init
  sigma <- ssm$G %*% ssm$Q %*% t(ssm$G)
  dim <- ssm$dim; Y <- ssm$Y
  
  # difference here 
  c <- ssm$c; delta <- ssm$delta
  
  pre_mu_1 <- solve(F %*% F + solve(sigma_init)%*%sigma) %*% F
  sigma_1 <- solve(F %*% (solve(sigma) %*% F) + solve(sigma_init))
  pre_mu_j <- solve(F %*% F + diag(1,dim(F)[1])) %*% F
  sigma_j <- solve(F %*% (solve(sigma) %*% F) + solve(sigma))
  pre_mu_T <- F; sigma_T <- sigma
  
  ##### MHMCMC autoregressive ######
  # N <- 1000; e <- 0.5 # ranges from -1 to 1 
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
      # hastings <- dmvnorm(Y[j,],H %*% x_j, R)/dmvnorm(Y[j,],H %*% X_mcmc[i-1,j,], R)
      lhastings <- 0
      for(p in 1:dim){
        num <- dpois(Y[j,p],exp(c[p]+delta[p]*x_j[p]),log=TRUE)
        denom <- dpois(Y[j,p],exp(c[p]+delta[p]*X_mcmc[i-1,j,p]),log=TRUE)
        lhastings <- lhastings + num - denom
      }
      alpha <- min(1,exp(lhastings))
      U <- runif(1,0,1)
      if(alpha>U){
        X_mcmc[i,j,] <- x_j
      } else {
        X_mcmc[i,j,] <- X_mcmc[i-1,j,]
      }
    }
  }
  close(pb)
  return(list(X_mcmc=X_mcmc,N=N,e=e))
}

mcmcGaussianSSM<- function(N,es,ssm,obs='Gaussian',init=NULL){
  require(MASS)
  # autoregressive proposal
  stopifnot((obs=='Gaussian')|(obs=='Poisson'))
  if(obs=='Gaussian'){
    require(Rfast)
    H <- ssm$H; R <- ssm$R; Q <- ssm$Q
  } else if(obs=='Poisson'){
    c <- ssm$c; delta <- ssm$delta
  }
  # extract ssm values
  F <- ssm$F; T <- ssm$T; Q <- ssm$Q
  mu_init <- ssm$mu_init; sigma_init <- ssm$sigma_init
  sigma <- ssm$G %*% ssm$Q %*% t(ssm$G)
  dim <- ssm$dim; Y <- ssm$Y
  
  # difference here 
  
  pre_mu_1 <- solve(F %*% F + solve(sigma_init)%*%sigma) %*% F
  sigma_1 <- solve(F %*% (solve(sigma) %*% F) + solve(sigma_init))
  pre_mu_j <- solve(F %*% F + diag(1,dim(F)[1])) %*% F
  sigma_j <- solve(F %*% (solve(sigma) %*% F) + solve(sigma))
  pre_mu_T <- F; sigma_T <- sigma
  
  ##### MHMCMC autoregressive ######
  # N <- 1000; e <- 0.5 # ranges from -1 to 1 
  X_mcmc <- array(0,dim=c(N,T,dim))
  if(is.null(init)){
    X_mcmc[1,,] <- mvrnorm(T,mu_init,sigma_init) 
  } else {
    X_mcmc[1,,] <- init
  }
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
      e <- es[i%%length(es)+1] # alternate
      x_j <- mu_j+sqrt(1-e^2)*(X_mcmc[i-1,j,]-mu_j)+e*L%*%z
      # transition probability

      if(obs=='Gaussian'){
        hastings <- dmvnorm(Y[j,],H %*% x_j, R)/dmvnorm(Y[j,],H %*% X_mcmc[i-1,j,], R)
      } else if(obs=='Poisson'){
        lhastings <- 0
        for(p in 1:dim){
          num <- dpois(Y[j,p],exp(c[p]+delta[p]*x_j[p]),log=TRUE)
          denom <- dpois(Y[j,p],exp(c[p]+delta[p]*X_mcmc[i-1,j,p]),log=TRUE)
          lhastings <- lhastings + num - denom
        }
        hastings = exp(lhastings)
      }

      alpha <- min(1,hastings)
      U <- runif(1,0,1)
      if(alpha>U){
        X_mcmc[i,j,] <- x_j
      } else {
        X_mcmc[i,j,] <- X_mcmc[i-1,j,]
      }
    }
  }
  close(pb)
  return(list(X_mcmc=X_mcmc,N=N,e=e,init=init))
}

