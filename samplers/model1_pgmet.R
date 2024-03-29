source('../../samplers/model1_pgbs.R') # require its forward and backward definitions

mcmc_step<- function(X_current,Y,N,es,T,dim,mus,sigmas_L,c,delta,thin.factor=10,seed=NULL){
  require(MASS)
  
  ##### MHMCMC autoregressive ######
  X_sample <- array(0,dim=c(N/thin.factor,T,dim))
  X_sample[1,,] <- X_current
  if(!is.null(seed)){
    set.seed(seed)
  }
  for(i in 2:N){
    U <- runif(T)
    for(j in 1:T){
      # find the mu 
      if(j==1) {
        # mu_j <- pre_mu_1 %*% X_sample[i-1,2,]
        mu_j <- mus[[1]] %*% X_current[2,]
        L <- sigmas_L[[1]]
      } else if(j==T) {
        mu_j <- mus[[3]] %*% X_current[T-1,]
        L <- sigmas_L[[3]]
      } else {
        # mu_j <- pre_mu_j %*% (X_sample[i,j-1,]+X_sample[i-1,j+1,])
        mu_j <- mus[[2]] %*% (X_current[j-1,]+X_current[j+1,])
        L <- sigmas_L[[2]]
      }
      
      # autoregressive update
      z <- rnorm(dim)
      e <- es[i%%length(es)+1] # alternate
      x_j <- mu_j+sqrt(1-e^2)*(X_current[j,]-mu_j)+e*L%*%z
      
      lnum <- sum(dpois(Y[j,],exp(c+delta*x_j),log=TRUE))
      # ldenom <- sum(dpois(Y[j,],exp(c+delta*X_sample[i-1,j,]),log=TRUE))
      ldenom <- sum(dpois(Y[j,],exp(c+delta*X_current[j,]),log=TRUE))
      lhastings <- lnum - ldenom
      hastings = exp(lhastings)
      
      alpha <- min(1,hastings)
      if(alpha>U[j]){
        X_current[j,] <- x_j
      } #else unchanged
    }
    
    if(i%%thin.factor==0){
      X_sample[i/thin.factor,,] <- X_current
    }
  }
  return(X_current) # return only the last sample
  # return(list(X_sample=X_sample,N=N,thin.factor=thin.factor,es=es,init=init))
}

pgmetModel1 <- function(ssm,N,L,es,N.mcmc=10,init=NULL,seed=NULL,return.weight=FALSE){
  require(MASS)
  #poisson observation and gaussian latent process 
  mu_init <- ssm$mu_init
  sigma_init <- ssm$sigma_init
  F <- ssm$F
  c <- ssm$c
  delta <- ssm$delta
  T <- ssm$T
  sigma <- ssm$sigma
  sigma_U <- ssm$sigma_U
  sigma_L <- ssm$sigma_L
  Y <- ssm$Y
  dim <- ssm$dim
  
  # set up for mcmc
  pre_mu_1 <- solve(F %*% F + solve(sigma_init)%*%sigma) %*% F
  sigma_1 <- solve(F %*% (solve(sigma) %*% F) + solve(sigma_init))
  pre_mu_j <- solve(F %*% F + diag(1,dim(F)[1])) %*% F
  sigma_j <- solve(F %*% (solve(sigma) %*% F) + solve(sigma))
  pre_mu_T <- F
  sigma_T <- sigma
  
  sigma_1_L <- t(chol(sigma_1))
  sigma_j_L <- t(chol(sigma_j))
  sigma_T_L <- t(chol(sigma_T))
  
  mus <- list(pre_mu_1,pre_mu_j,pre_mu_T)
  sigmas_L <- list(sigma_1_L,sigma_j_L,sigma_T_L)
  
  if(is.null(ssm$sigma_inv)){
    print('old ssm object, computing inverse for sigma')
    sigma_inv <- chol2inv(sigma_U)
  } else {
    sigma_inv <- ssm$sigma_inv
  }

  # intialise sample matrix
  X_sample <- array(0,dim=c(4*N+1,T,dim))
  if(return.weight){
    W <- array(logical(0),dim=c(2*N,T,L))
    lW <- array(logical(0),dim=c(2*N,T,L))
  }
  if(is.null(init)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    X_sample[1,,] <- mvrnorm(T,mu_init,sigma_init)
  } else {
    X_sample[1,,] <- init 
  }
  X_current <- X_sample[1,,]
  
  pb <- txtProgressBar(min=0,max=N*4,title="pgbs",style=3)
  for(n in seq(2,4*N,4)){
    ### Step1: pgbs with forward sequence ###
    forward_results <- forward_step(X_current,Y,T,L,dim,F,sigma_U,mu_init,sigma_init,c,delta)
    if(return.weight){
      W[n/2,,] <- forward_results$W
      lW[n/2,,] <- forward_results$lW
    }
    X_sample[n,,] <- backward_step(forward_results,T,L,F,sigma_inv)
    
    ### Step2: mcmc 10 steps ###
    X_sample[n+1,,] <- mcmc_step(X_sample[n,,],Y,N.mcmc,es,T,dim,mus,sigmas_L,c,delta,thin.factor=1)
    
    ### Step3: pgbs reversed sequence ###
    forward_results <- forward_step(X_sample[n+1,seq(T,1,-1),],Y[seq(T,1,-1),],T,L,dim,F,sigma_U,mu_init,sigma_init,c,delta)
    if(return.weight){
      W[n/2+1,seq(T,1,-1),] <- forward_results$W
      lW[n/2+1,seq(T,1,-1),] <- forward_results$lW
    }
    X_sample[n+2,seq(T,1,-1),] <- backward_step(forward_results,T,L,F,sigma_inv) # reverse the index when setting
    
    ### Step4: mcmc 10 steps ###
    X_sample[n+3,,] <- mcmc_step(X_sample[n+2,,],Y,N.mcmc,es,T,dim,mus,sigmas_L,c,delta,thin.factor=1)
    
    ### update current and terminate ###
    X_current <- X_sample[n+3,,] # forward sequence, updated once for next iteration only
    setTxtProgressBar(pb, n)
  }
  setTxtProgressBar(pb, 4*N)
  close(pb)
  if(return.weight){
    return(list(X_sample = X_sample[-1,,],N=N,init=init,W=W,lW=lW,seed=seed))
  }
  return(list(X_sample = X_sample[-1,,],N=N,init=init,seed=seed))
}