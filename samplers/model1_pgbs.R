ld_pois_model1 <- function(x,y,c,delta){
  # note: independent components, specific to model 1
  ldense <- 0 
  lambda <- c+delta*x
  ldense <- dpois(y,exp(lambda),log=TRUE)
  return(sum(ldense))
}

r_model1_transition <- function(x_from,F,sigma_L){
  return(F%*%x_from+sigma_L%*%rnorm(length(x_from)))
}

# d_model1_transition <- function(x_from,x_to,F,sigma_U){
#   d <- length(x_to)
#   det_U <- sum(diag(sigma_U))
#   # det_U <- tr(sigma_U)
#   inv_sigma <- chol2inv(sigma_U)
#   diff <- x_to-F%*%x_from
#   dense <- (2*pi)^(-d/2) / det_U * exp(-1/2*t(diff)%*%inv_sigma%*%diff)
#   return(dense)
# }

ld_model1_transition <- function(x_from,x_to,F,sigma_U){
  inv_sigma <- chol2inv(sigma_U)
  diff <- x_to-F%*%x_from
  return(-1/2*t(diff)%*%inv_sigma%*%diff)
}


pgbsModel1 <- function(ssm,N,L,init){
  #poisson observation and gaussian latent process 
  X_sample <- array(0,dim=c(N+1,T,dim))
  X_sample[1,,] <- init 
  mu_init <- ssm$mu_init; sigma_init <- ssm$sigma_init
  F <- ssm$F; sigma <- ssm$sigma; c <- ssm$c; delta <- ssm$delta
  T <- ssm$T
  sigma_U <- ssm$sigma_U; sigma_L <- ssm$sigma_L
  pb <- txtProgressBar(min=0,max=N,title="pgbs",style=3)
  for(n in 2:N){
    setTxtProgressBar(pb, n)
    # print(paste("progress: ", n, '/', N,sep=''))
    x_pool <- array(0,dim=c(T,L,dim))
    W <- matrix(0,nrow=T,ncol=L)
    
    ### Forward pass ### 
    ### t=1 ###
    x_pool[1,1,] <- X_sample[n-1,1,]
    x_pool[1,2:L,] <- rmvnorm(L-1,mu_init,sigma_init)
    W[1,] <- apply(x_pool[1,,],MARGIN=1,FUN=mvdpois_model1,y=ssm$Y[1,],c=c,delta=delta)
    W[1,] <- W[1,]/sum(W[1,])
    
    ### t=2,3,.. ###
    for(t in 2:T){
      x_pool[t,1,] <- X_sample[n-1,t,]
      anc <- sample(1:L,L-1,replace=TRUE,prob=W[t-1,])
      x_pool[t,2:L,] <- apply(x_pool[t-1,anc,],MARGIN=1,FUN=r_model1_transition,F=F,sigma_L=sigma_L)
      W[t,] <- apply(x_pool[t,,],MARGIN=1,FUN=ld_pois_model1,y=ssm$Y[t,],c=c,delta=delta)
      W[t,] <- exp(W[t,]-max(W[t,]))/sum(exp(W[t,]-W[t,]))
    }
    
    ### backward pass ###
    l_T <- sample(1:L,1,replace=TRUE,prob=W[T,])
    X_sample[n,T,] <- x_pool[T,l_T,]
    for(t in (T-1):1){
      prob <- apply(x_pool[t,,],MARGIN=1,FUN=ld_model1_transition,x_to=X_sample[n,t+1,],F=F,sigma_U=sigma_U)
      prob <- exp(prob) * W[t,]
      prob <- prob/sum(prob)
      # setting new x
      X_sample[n,t,] <- x_pool[t,sample(1:L,1,prob=prob),]
    }
  }
  return(X_sample[2:N,,])
  # modify to include the other sampling parameters 
}

pgbsModel1_batch <- function(ssm,N,L,init,batch){
  #poisson observation and gaussian latent process 
  X_sample <- array(0,dim=c(N*batch+1,T,dim))
  X_sample[1,,] <- init 
  mu_init <- ssm$mu_init; sigma_init <- ssm$sigma_init
  F <- ssm$F; sigma <- ssm$sigma; c <- ssm$c; delta <- ssm$delta
  T <- ssm$T
  sigma_U <- ssm$sigma_U; sigma_L <- ssm$sigma_L
  pb <- txtProgressBar(min=0,max=batch*N,title="pgbs",style=3)
  for(n in seq(2,(batch*N),batch)){
    
    # print(paste("progress: ", n, '/', N,sep=''))
    x_pool <- array(0,dim=c(T,L,dim))
    W <- matrix(0,nrow=T,ncol=L)
    
    ### Forward pass ### 
    ### t=1 ###
    x_pool[1,1,] <- X_sample[n-1,1,]
    x_pool[1,2:L,] <- rmvnorm(L-1,mu_init,sigma_init)
    W[1,] <- apply(x_pool[1,,],MARGIN=1,FUN=ld_pois_model1,y=ssm$Y[1,],c=c,delta=delta)
    W[1,] <- exp(W[1,]-max(W[1,]))/sum(exp(W[1,]-W[1,]))
    
    ### t=2,3,.. ###
    for(t in 2:T){
      x_pool[t,1,] <- X_sample[n-1,t,]
      anc <- sample(1:L,L-1,replace=TRUE,prob=W[t-1,])
      x_pool[t,2:L,] <- apply(x_pool[t-1,anc,],MARGIN=1,FUN=r_model1_transition,F=F,sigma_L=sigma_L)
      W[t,] <- apply(x_pool[t,,],MARGIN=1,FUN=ld_pois_model1,y=ssm$Y[t,],c=c,delta=delta)
      W[t,] <- exp(W[t,]-max(W[t,]))/sum(exp(W[t,]-W[t,]))
    }
    
    ### backward pass ###
    l_T <- sample(1:L,batch,replace=TRUE,prob=W[T,])
    X_sample[n:(n+batch-1),T,] <- x_pool[T,l_T,]
    for(t in (T-1):1){
      prob <- apply(x_pool[t,,],MARGIN=1,FUN=ld_model1_transition,x_to=X_sample[n,t+1,],F=F,sigma_U=sigma_U)
      prob <- exp(prob) * W[t,]
      prob <- prob/sum(prob)
      # setting new x
      X_sample[n:(n+batch-1),t,] <- x_pool[t,sample(1:L,batch,replace=TRUE,prob=prob),]
      # print(x_pool[t,sample(1:L,2,replace=TRUE,prob=prob),])
    }
    setTxtProgressBar(pb, n)
  }
  close(pb)
  return(list(X_sample = X_sample[-1,,],N=N,init=init))
  # modify to include the other sampling parameters 
}
