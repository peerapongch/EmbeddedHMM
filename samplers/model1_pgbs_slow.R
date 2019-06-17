source('./model1_pgbs.R')
pgbsModel1_slow <- function(ssm,N,L,x_init){
  #poisson observation and gaussian latent process 
  x_new <- array(0,dim=c(N,T,dim))
  mu_init <- ssm$mu_init; sigma_init <- ssm$sigma_init
  F <- ssm$F; sigma <- ssm$sigma; c <- ssm$c; delta <- ssm$delta
  T <- ssm$T
  for(n in 1:N){
    # print(paste("progress: ", n, '/', N,sep=''))
    x_pool <- array(0,dim=c(T,L,dim))
    W <- matrix(0,nrow=T,ncol=L)
    
    ### Forward pass ### 
    ### t=1 ###
    x_pool[1,1,] <- x_init[1,]
    x_pool[1,2:L,] <- rmvnorm(L-1,mu_init,sigma_init)
    
    for(l in 1:L){
      W[1,l] <- mvdpois_model1(x_pool[1,l,],ssm_poisson$Y[1,],c,delta)
    }
    
    W[1,] <- W[1,]/sum(W[1,])
    
    ### t=2,3,.. ###
    for(t in 2:T){
      x_pool[t,1,] <- x_init[t,]
      W[t,1] <- mvdpois_model1(ssm$Y[t,],x_pool[t,1,],c,delta)
      
      anc <- sample(1:L,L,replace=TRUE,prob=W[t-1,])
      for(l in 2:L){
        x_pool[t,l,] <- rmvnorm(1,F %*% x_pool[t-1,anc[l],],sigma)
        W[t,l] <- mvdpois_model1(x_pool[t,l,],ssm_poisson$Y[t,],c,delta)
      }
      W[t,] <- W[t,]/sum(W[t,])
    }
    
    ### backward pass ###
    # x_new <- matrix(0,nrow=T,ncol=dim)
    l_T <- sample(1:L,1,replace=TRUE,prob=W[T,])
    x_new[n,T,] <- x_pool[T,l_T,]
    for(t in (T-1):1){
      prob <- rep(0,L)
      for(l in 1:L){
        prob[l] <- dmvnorm(x_new[n,t+1,],F%*%x_pool[t,l,],sigma)
      }
      prob <- prob * W[t,]; prob <- prob/sum(prob)
      # setting new x
      x_new[n,t,] <- x_pool[t,sample(1:L,1,prob=prob),]
    }
  }
  return(x_new)
  # modify to include the other sampling parameters 
}
