mvdpois_model1 <- function(x,y,c,delta){
  # note: independent components, specific to model 1
  ldense <- 0 
  for(i in 1:length(y)){
    ldense <- ldense + dpois(y[i],exp(c[i]+delta[i]*x[i]),log=TRUE)
  }
  return(exp(ldense))
}

r_model1_transition <- function(x_from,F,sigma_L){
  return(F%*%x_from+sigma_L%*%rnorm(length(x_from)))
}

d_model1_transition <- function(x_from,x_to,F,sigma_U){
  d <- length(x_to)
  det_U <- sum(diag(sigma_U))
  # det_U <- tr(sigma_U)
  inv_sigma <- chol2inv(sigma_U)
  diff <- x_to-F%*%x_from
  dense <- (2*pi)^(-d/2) / det_U * exp(-1/2*t(diff)%*%inv_sigma%*%diff)
  return(dense)
}


pgbsModel1 <- function(ssm,N,L,x_init){
  #poisson observation and gaussian latent process 
  x_new <- array(0,dim=c(N,T,dim))
  mu_init <- ssm$mu_init; sigma_init <- ssm$sigma_init
  F <- ssm$F; sigma <- ssm$sigma; c <- ssm$c; delta <- ssm$delta
  T <- ssm$T
  sigma_U <- chol(sigma)
  sigma_L <- t(sigma_U)
  for(n in 1:N){
    # print(paste("progress: ", n, '/', N,sep=''))
    x_pool <- array(0,dim=c(T,L,dim))
    W <- matrix(0,nrow=T,ncol=L)
    
    ### Forward pass ### 
    ### t=1 ###
    x_pool[1,1,] <- x_init[1,]
    x_pool[1,2:L,] <- rmvnorm(L-1,mu_init,sigma_init)
    W[1,] <- apply(x_pool[1,,],MARGIN=1,FUN=mvdpois_model1,y=ssm$Y[1,],c=c,delta=delta)
    W[1,] <- W[1,]/sum(W[1,])
    
    ### t=2,3,.. ###
    for(t in 2:T){
      x_pool[t,1,] <- x_init[t,]
      anc <- sample(1:L,L-1,replace=TRUE,prob=W[t-1,])
      x_pool[t,2:L,] <- apply(x_pool[t-1,anc,],MARGIN=1,FUN=r_model1_transition,F=F,sigma_L=sigma_L)
      W[t,] <- apply(x_pool[t,,],MARGIN=1,FUN=mvdpois_model1,y=ssm$Y[t,],c=c,delta=delta)
      W[t,] <- W[t,]/sum(W[t,])
    }
    
    ### backward pass ###
    # x_new <- matrix(0,nrow=T,ncol=dim)
    l_T <- sample(1:L,1,replace=TRUE,prob=W[T,])
    x_new[n,T,] <- x_pool[T,l_T,]
    for(t in (T-1):1){
      prob <- rep(0,L)
      
      # optimisation: pre-compute the difference
      prob <- apply(x_pool[t,,],MARGIN=1,FUN=d_model1_transition,x_to=x_new[n,t+1,],F=F,sigma_U=sigma_U)
      
      # for(l in 1:L){
      #   prob[l] <- dmvnorm(x_new[n,t+1,],F%*%x_pool[t,l,],sigma)
      # }
      prob <- prob * W[t,]; prob <- prob/sum(prob)
      # setting new x
      x_new[n,t,] <- x_pool[t,sample(1:L,1,prob=prob),]
    }
  }
  return(x_new)
  # modify to include the other sampling parameters 
}
