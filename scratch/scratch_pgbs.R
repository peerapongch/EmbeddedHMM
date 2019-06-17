load('./data/ssm_model1_poisson_1.RData')
load('./data/ssm_model1_poisson_1_env.RData')
library(mvtnorm)
mvdpois_model1 <- function(y,x,c,delta){
  # note: independent components, specific to model 1
  ldense <- 0 
  for(i in 1:length(y)){
    ldense <- ldense + dpois(y[i],exp(c[i]+delta[i]*x[i]),log=TRUE)
  }
  return(exp(ldense))
}

L <- 100
x_init <- matrix(0,nrow=T,ncol=dim) # init ofc refers to first iteration, replace this with x_current for later iterations



x_pool <- array(0,dim=c(T,L,dim))
W <- matrix(0,nrow=T,ncol=L)

### t=1 ###
x_pool[1,1,] <- x_init[1,]
x_pool[1,2:L,] <- mvrnorm(L-1,mu_init,sigma_init)
for(l in 1:L){
  W[1,l] <- mvdpois_model1(ssm_poisson$Y[1,],x_pool[1,l,],c,delta)
}
W[1,] <- W[1,]/sum(W[1,])

### t=2,3,.. ###
for(t in 2:T){
  x_pool[t,1,] <- x_init[t,]
  W[t,1] <- mvdpois_model1(ssm_poisson$Y[t,],x_pool[t,1,],c,delta)
  anc <- sample(1:L,L,replace=TRUE,prob=W[1,])
  for(l in 2:L){
    x_pool[t,l,] <- mvrnorm(1,F %*% x_pool[t-1,anc[l],],sigma)
    W[t,l] <- mvdpois_model1(ssm_poisson$Y[t,],x_pool[t,l,],c,delta)
  }
  W[t,] <- W[t,]/sum(W[t,])
}

### backward pass ###
x_new <- matrix(0,nrow=T,ncol=dim)
l_T <- sample(1:L,1,replace=TRUE,prob=W[T,])
x_new[T,] <- x_pool[T,l_T,]
for(t in (T-1):1){
  prob <- rep(0,L)
  for(l in 1:L){
    prob[l] <- dmvnorm(x_new[t+1,],F%*%x_pool[t,l,],sigma)
  }
  prob <- prob * W[t,]; prob <- prob/sum(prob)
  # setting new x
  x_new[t,] <- x_pool[t,sample(1:L,1,prob=prob),]
}
