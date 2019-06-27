autoregressive_update <- function(this_X_current,mu,Y,t,l,l_current,F,c,delta,sq_term,Us,last){
  # autoregressive_update <- function(this_X_current,mu,Y,t,l,l_current,F,sigma_L,c,delta,es,es_2,zs,Us){
  # print(F)
  # print(X_pool[t,l_current,])
  # mu <- F %*% X_pool[t,l_current,]
  mu_current <- mu[l_current,]
  # X_new <- mu_current + sqrt(1-es_2[l])*(this_X_current-mu_current)+es[l]*sigma_L%*%zs[l,]
  X_new <- mu_current + sq_term[l]*(this_X_current-mu_current)+last[l,]
  # hastings ratio
  lnum <- sum(dpois(Y[t,],exp(c+delta*X_new),log=TRUE))
  ldenom <- sum(dpois(Y[t,],exp(c+delta*this_X_current),log=TRUE))
  lhastings <- lnum - ldenom
  hastings <- exp(lhastings)
  alpha <- min(1,hastings)
  if(alpha>Us[1,l]){
    return(X_new)
  }
  return(this_X_current)
}

shift_update <- function(this_X_current,X_pool,Y,t,l,l_current,F,sigma_L,c,delta,Us,l_new){
  # shift update: then update x and l using shift update 
  # update of l UAR(1:L) done above since independent of l
  X_new <- this_X_current + F %*% (X_pool[t-1,l_new[l],]-X_pool[t-1,l_current,])
  lnum <- sum(dpois(Y[t,],exp(c+delta*X_new),log=TRUE))
  ldenom <- sum(dpois(Y[t,],exp(c+delta*this_X_current),log=TRUE))
  lhastings <- lnum - ldenom
  hastings <- exp(lhastings)
  alpha <- min(1,hastings)
  if(alpha>Us[2,l]){
    return(list(X_new=X_new,l_new=l_new[l]))
  } #else unchanged
  return(list(X_new=this_X_current,l_new=l_current))
}

forward_pool <- function(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,F,sigma_inv,c,delta,reversed){
  # metropolis
  # sequential update
  X_pool <- array(0,dim=c(T,L,dim))
  
  # start
  es <- runif(L,0.1,0.4) # hardcoded
  es_2 <- es^2
  zs <- matrix(rnorm(L*dim),ncol=dim,nrow=L)
  Us <- runif(L)
  for(l in 1:L){
    # autoregressive update 
    X_new <- mu_init + sqrt(1-es_2[l])*(X_current[1,]-mu_init) + es[l]*sigma_init_L%*%zs[l,]
    # hastings ratio
    lnum <- sum(dpois(Y[1,],exp(c+delta*X_new),log=TRUE))
    ldenom <- sum(dpois(Y[1,],exp(c+delta*X_current[1,]),log=TRUE))
    lhastings <- lnum - ldenom
    hastings <- exp(lhastings)
    
    alpha <- min(1,hastings)
    if(alpha>Us[l]){
      X_current[1,] <- X_new
    } #else unchanged
    
    X_pool[1,l,] <- X_current[1,]
  }
  
  # then for t>1
  for(t in 2:T){
    # l needs to be sampled in sequence :( 
    diff <- X_current[t,] - t(X_pool[t-1,,] %*% F) # for symmetric F
    l_lprob <- diag(-1/2*t(diff)%*%sigma_inv%*%diff)
    prob <- exp(l_lprob-max(l_lprob)); prob <- prob/sum(prob)
    l_current <- sample(1:L,1,prob=prob)
    
    es <- runif(L,0.1,0.4) # hardcoded
    es_2 <- es^2
    zs <- matrix(rnorm(L*dim),ncol=dim,nrow=L) # note dim reversed from above
    Us <- matrix(runif(L*2),nrow=2,ncol=L)
    l_new <- sample(1:L,L,replace=TRUE)
    
    # experiment to cut cost of setting current to pool: worked
    this_X_current <- X_current[t,]
    this_X_pool <- matrix(logical(0),nrow=L,ncol=dim)
    mu <- X_pool[t-1,,] %*% F
    
    sq_term <- sqrt(1-es^2)
    last <- es*zs%*%sigma_L
    for(l in 1:L){
      if(!reversed){
        # AR then shift
        X_new <- autoregressive_update(this_X_current,mu,Y,t,l,l_current,F,c,delta,sq_term,Us,last)
        # count and then 
        this_X_current <- X_new
        
        shift_out <- shift_update(this_X_current,X_pool,Y,t,l,l_current,F,sigma_L,c,delta,Us,l_new)
        # count and then update
        this_X_current <- shift_out$X_new
        l_current <- shift_out$l_new
        # set 
        this_X_pool[l,] <- this_X_current
        # X_pool[t,l,] <- this_X_current  
      } else {
        # shift then AR
        shift_out <- shift_update(this_X_current,X_pool,Y,t,l,l_current,F,sigma_L,c,delta,Us,l_new)
        # count and then update
        this_X_current <- shift_out$X_new
        l_current <- shift_out$l_new
        
        # X_new <- autoregressive_update(this_X_current,mu,Y,t,l,l_current,F,sigma_L,c,delta,es,es_2,zs,Us)
        X_new <- autoregressive_update(this_X_current,mu,Y,t,l,l_current,F,c,delta,sq_term,Us,last)
        # count and then 
        this_X_current <- X_new
        
        # set 
        this_X_pool[l,] <- this_X_current
        # X_pool[t,l,] <- this_X_current
      }
    }
    X_pool[t,,] <- this_X_pool
    X_current[t,] <- this_X_current
  } 
  return(X_pool)
}

backward_sampling <- function(X_pool,L,T,dim,F,sigma_inv){
  # backward sampling
  X_new <- matrix(0,nrow=T,ncol=dim)
  X_new[T,] <- X_pool[T,sample(1:L,1),] # since UAR on 1 to L @T 
  for(t in (T-1):1){
    # batch calculate the transition probability to the chosen
    diff <- X_new[t+1,] - t(X_pool[t,,] %*% F) # also only for symmetric F 
    lprob <- diag(-1/2*t(diff)%*%sigma_inv%*%diff)
    prob <- exp(lprob-max(lprob)); prob <- prob/sum(prob)
    X_new[t,] <- X_pool[t,sample(1:L,1,prob=prob),]
  }
  return(X_new)
}

ehmmModel1 <- function(ssm,N,L,init=NULL,seed=NULL){
  # sampling epsilon values instead
  # setup
  require(MASS)
  #poisson observation and gaussian latent process 
  mu_init <- ssm$mu_init
  sigma_init <- ssm$sigma_init
  F <- ssm$F
  c <- ssm$c
  delta <- ssm$delta
  T <- ssm$T
  sigma_U <- ssm$sigma_U
  sigma_L <- ssm$sigma_L
  Y <- ssm$Y
  dim <- ssm$dim
  # for backward compatibility 
  if(is.null(ssm$sigma_init_L)){
    print('old ssm object, computing Choleskey for sigma_init')
    sigma_init_L <- t(chol(sigma_init))
  } else {
    sigma_init_L <- ssm$sigma_init_L
  }
  
  if(is.null(ssm$sigma_inv)){
    print('old ssm object, computing inverse for sigma')
    sigma_inv <- chol2inv(sigma_U)
  } else {
    sigma_inv <- ssm$sigma_inv
  }
  
  # X_sample <- array(0,dim=c(2*N+1,T,dim))
  X_sample <- array(0,dim=c(2*N+1,T,dim))
  if(is.null(init)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    X_sample[1,,] <- mvrnorm(T,mu_init,sigma_init)
  } else {
    X_sample[1,,] <- init 
  }
  X_current <- X_sample[1,,]
  pb <- txtProgressBar(min=0,max=2*N,title="pgbs",style=3)
  for(i in seq(2,2*N,2)){
    # forward sequence
    X_pool <- forward_pool(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,F,sigma_inv,c,delta,reversed=FALSE)
    X_current <- backward_sampling(X_pool,L,T,dim,F,sigma_inv)
    X_sample[i,,] <- X_current
    
    # reversed sequence
    ## todo (modify forward_pool to detect reversed chain), and test forward
    X_pool <- forward_pool(X_current[seq(T,1,-1),],Y[seq(T,1,-1),],T,L,dim,mu_init,sigma_init_L,sigma_L,F,sigma_inv,c,delta,reversed=TRUE)
    X_current[seq(T,1,-1),] <- backward_sampling(X_pool,L,T,dim,F,sigma_inv)
    X_sample[i+1,,] <- X_current
    
    setTxtProgressBar(pb, i)
  }
  return(list(X_sample=X_sample[-1,,],N=N,L=L,init=init,seed=seed,X_pool=X_pool))
}