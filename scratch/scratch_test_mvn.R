x <- mvrnorm(1,mu_init,sigma_init)

system.time(
  replicate(10000,mvrnorm(1,F %*% mu_init,sigma_init))
)

system.time(
  replicate(10000,Rfast::rmvnorm(1,F %*% mu_init,sigma_init))
)


system.time(
  replicate(10000,Rfast::dmvnorm(x_from,mu_init,sigma_init))
)

system.time(
  replicate(10000,mvtnorm::dmvnorm(x_from,mu_init,sigma_init))
)

system.time(
  replicate(10000,mvnfast::dmvn(x_from,mu_init,sigma_init))
)

apply(x_pool[1,,],MARGIN=1,mvdpois_model1,y=ssm_poisson$Y[1,],c=c,delta=delta)
dim(x_pool[1,,])

### custom functions ###
x_from <- mvrnorm(1,mu_init,sigma_init)
custom_rmvnorm <- function(F,x_from,G){
  return(F%*%x_from+G%*%rnorm(length(x_from)))
}
system.time(
  replicate(10000,custom_rmvnorm(F,x_from,G))
)


custom_dmvnorm <- function(x,mu,sigma_chol){
  d <- length(x)
  det_chol <- sum(diag(sigma_chol))
  inv_chol <- solve(sigma_chol)
  inv_sigma <- t(inv_chol) %*% inv_chol
  diff <- x-mu
  dense <- (2*pi)^(-d/2) / det_chol * exp(-1/2*t(diff)%*%inv_sigma%*%diff)
  return(dense)
}

L <- t(chol(sigma_init))
system.time(
  replicate(10000,custom_dmvnorm(x_from,mu_init,L))
)




### density calculation speed up please!! ###
d_model1_transition <- function(x_from,x_to,F,sigma){
  return(dmvnorm(x_to,F %*% x_from,sigma))
}

d_model1_transition <- function(x,mu,sigma_chol){
  d <- length(x)
  det_chol <- sum(diag(sigma_chol))
  inv_chol <- solve(sigma_chol)
  inv_sigma <- t(inv_chol) %*% inv_chol
  diff <- x-mu
  dense <- (2*pi)^(-d/2) / det_chol * exp(-1/2*t(diff)%*%inv_sigma%*%diff)
  return(dense)
}

####
system.time(replicate(1000000,t(G)))
system.time(replicate(1000000,Rfast::transpose(G)))

system.time(replicate(100000,t(G)))
system.time(replicate(100000,t(x_init)))
