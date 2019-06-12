# implements stationary Gaussian SSM with stationary Gaussian latent process and stationary Gaussian observation process 

generateLatent_Gaussian <- function(T,dim,F,G,Q){
  
}

generateObservation_Gaussian <- function(T,X,H,R){
  # X being the latent states 
  # ambiguous how to make use of the dimensions 
  
}

generateGaussianGaussianSSM <- function(T,dim,mu_init,sigma_init,F,G,Q,H,R){
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
  # generate observations 
  Y <- matrix(0,ncol=dim,nrow=T)
  for(t in 1:T){
    Y[t,] <- mvrnorm(1, H %*% X[t,], R)
  }
  
  return(list(dim=dim,T=T,x0=x0,X=X,Y=Y,mu_init=mu_init,sigma_init=sigma_init,F=F,G=G,Q=Q,H=H,R=R))
}