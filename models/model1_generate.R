# implements stationary Gaussian SSM with stationary Gaussian latent process and stationary Gaussian observation process 

generateLatent_Gaussian <- function(T,dim,mu_init,sigma_init,F,G,Q){
  require(MASS)
  # init
  sigma <- G %*% Q %*% t(G)
  # begin
  x0 <- mvrnorm(1,mu_init,sigma_init)
  X <- matrix(0,ncol=dim,nrow=T)
  for(t in 1:T){
    if(t==1){
      X[t,] <- mvrnorm(1, F %*% x0, sigma)  
    } else {
      X[t,] <- mvrnorm(1, F %*% X[t-1,], sigma) 
    }
  }
  return(list(x0=x0,X=X))
}

generateObservation_Poisson <- function(T,dim,X,c,delta){
  # as in model 1 in the paper 
  Y <- matrix(0,ncol=dim,nrow=T)
  for(t in 1:T){
    for(j in 1:dim){
      Y[t,j] <- rpois(1,exp(c[j]+delta[j]*X[t,j]))
    }
  }
  return(Y)
}

generatePoissonGaussianSSM <- function(T,dim,mu_init,sigma_init,F,G,Q,c,delta){
  # call the above two functions 
  X <- generateLatent_Gaussian(T,dim,mu_init,sigma_init,F,G,Q)
  Y <- generateObservation_Poisson(T,dim,X$X,c,delta)
  sigma <- G %*% Q %*% t(G)
  sigma_U <- chol(sigma)
  sigma_L <- t(sigma_U)
  return(list(dim=dim,T=T,x0=X$x0,X=X$X,Y=Y,mu_init=mu_init,sigma_init=sigma_init,
              F=F,G=G,Q=Q,c=c,delta=delta,sigma=sigma,sigma_U=sigma_U,sigma_L=sigma_L))
}

generateObservation_Gaussian <- function(T,X,H,R){
  # X being the latent states 
  # ambiguous how to make use of the dimensions 
  print('not implemented')
}

generateGaussianGaussianSSM <- function(T,dim,mu_init,sigma_init,F,G,Q,H,R){
  require(MASS)
  # init
  sigma <- G %*% Q %*% t(G)
  sigma_U <- chol(sigma)
  sigma_L <- t(sigma_U)
  # begin
  x0 <- mvrnorm(1,mu_init,sigma_init)
  X <- matrix(0,ncol=dim,nrow=T)
  for(t in 1:T){
    if(t==1){
      X[t,] <- mvrnorm(1, F %*% x0, sigma)  
    } else {
      X[t,] <- mvrnorm(1, F %*% X[t-1,], sigma) 
    }
  }
  # generate observations 
  Y <- matrix(0,ncol=dim,nrow=T)
  for(t in 1:T){
    Y[t,] <- mvrnorm(1, H %*% X[t,], R)
  }
  
  return(list(dim=dim,T=T,x0=x0,X=X,Y=Y,mu_init=mu_init,sigma_init=sigma_init,
              F=F,G=G,Q=Q,H=H,R=R,sigma=sigma,sigma_U=sigma_U,sigma_L=sigma_L))
}