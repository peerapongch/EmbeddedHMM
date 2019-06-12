# model 1 as specified in the paper 
# but with gaussian observation model instaed of poisson 

# settings for the hidden model
dim_x <- 3
rho <- 0.7 
phi <- 0.9
phis <- rep(phi,dim_x)
#phis <- runif(3,0,1)
F <- diag(phis)


sigma <- matrix(rho,ncol=dim_x,nrow=dim_x)
for(i in 1:dim_x){
  sigma[i,i] <- 1
}

# Miscelleneous
G <- t(chol(sigma))
R <- diag(1,dim_x)
Q <- diag(1,dim_x)

mu_init <- rep(0,dim_x)
v <- 1/sqrt((1-phis^2))
sigma_init <- v %*% t(v)
for(i in 1:dim_x){
  for(j in 1:dim_x){
    if(i!=j){
      sigma_init[i,j] <- sigma_init[i,j] * rho
    }
  }
}

# settings for the observation model 
delta <- 0.6
dim_y <- dim_x
deltas <- rep(delta,dim_y)
H <- diag(deltas)
sd_y <- 1
sigma_y <- diag(sd_y,dim_y,dim_y)

# sample: latent
library(MASS)
T <- 250
x1 <- mvrnorm(1,mu_init,sigma_init)
X <- matrix(0,ncol=dim_x,nrow=n)
for(t in 1:T){
  if(t==1){
    X[t,] <- mvrnorm(1, F %*% x1, sigma)  
  } else {
    X[t,] <- mvrnorm(1, F %*% X[t-1,], sigma) 
  }
}
plot(X[,1],X[,2])

# sample: observation
deltas <- rep(0.6,dim_y)
H <- diag(deltas)
Y <- matrix(0,ncol=dim_y,nrow=n)
for(t in 1:n){
  Y[t,] <- mvrnorm(1, H %*% X[t,], sigma_y)
}

# comparison
par(mfrow=c(2,1))
plot(Y[,1])
plot(X[,1])
par(mfrow=c(1,1))
