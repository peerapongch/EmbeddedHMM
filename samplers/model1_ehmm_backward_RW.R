rw_update_x_oat_sys <- function(this_X_current,X_pool,Y,t,l,l_current,rw_const,sigma_inv,Us){
  # one at a time systematic scan randomwalk mhmcmc 
  # specificially for backward pool state sampling
  dim <- length(this_X_current)
  X_current <- this_X_current

  diff_current <- X_pool[t+1,l_current,] - F %*% X_current
  ldenom <- -1/2*t(diff_current)%*%sigma_inv%*%diff_current

  proposed <- rw_const*(2*runif(dim)-1)+X_current

  for(j in 1:dim){
    X_new <- X_current 
    X_new[j] <- proposed[j]

    diff_new <-  X_pool[t+1,l_current,] - F %*% X_new
    lnum <- -1/2*t(diff_new)%*%sigma_inv%*%diff_new
    
    hastings <- exp(lnum - ldenom)
    alpha <- min(1,hastings)
    # print(alpha)
    # print(hastings)
    # print(Us[1,l])
    # print('===')
    tryCatch({
      if(alpha>Us[1,l]){ 
        X_current <- X_new
        ldenom <- lnum
      }
    }, error = function(err){
      print('------------------')
      print(proposed)
      print('==')
      
      print(X_new)
      print(lnum)
      print('------------------')
    }
    )
  }
  return(X_current)
}

# normal_update <- function(this_X_current,X_pool,Y,t,l,l_current,rw_const,c,delta,sigma_inv,Us,last){
#   # this is autoregressive update without the nice cancellation in the hasting ratio 
#   dim <- length(this_X_current)
#   X_new <- this_X_current + last[l,]

#   Xs <- matrix(logical(0),nrow=dim,ncol=2)
#   Xs[,1] <- this_X_current
#   Xs[,2] <- X_new

#   diff <- as.vector(X_pool[t+1,l_current,]) - F %*% Xs
#   lprob <- diag(-1/2 * t(diff) %*% sigma_inv %*% diff)
#   lnum  <-  lprob[1]
#   ldenom <- lprob[2]

#   hastings <- exp(lnum-ldenom)
#   alpha <- min(1,hastings)
#   if(alpha>Us[1,l]){
#     return(X_new)
#   }
#   return(this_X_current)
# }

uar_update_l <- function(this_X_current,dim_index,X_pool,Y,t,l,l_current,c,delta,sigma_inv,l_new,Us){
  # propose l UAR then accept reject mh style
  lprob_obs_new <- sum(dpois(Y[t,],exp(c+delta*X_pool[t+1,l_new[l],]),log=TRUE))
  lprob_obs_current <- sum(dpois(Y[t,],exp(c+delta*X_pool[t+1,l_current,]),log=TRUE))

  # diff_new <- X_pool[t+1,l_new[l],] - F %*% this_X_current
  # diff_current <- X_pool[t+1,l_current,] - F %*% this_X_current

  diff <- t(X_pool[t+1,c(l_current,l_new[l]),]) - as.vector(F %*% this_X_current)
  lprob_trans <- diag(-1/2*t(diff)%*%sigma_inv%*%diff)

  lprob_trans_current <- lprob_trans[1]
  lprob_trans_new <- lprob_trans[2]
  # lprob_trans_new <- -1/2*t(diff_new)%*%sigma_inv%*%diff_new
  # lprob_trans_current <- -1/2*t(diff_current)%*%sigma_inv%*%diff_current

  hastings <- exp(lprob_trans_new+lprob_obs_new-lprob_trans_current-lprob_obs_current)
  alpha <- min(1,hastings)
  if(alpha>Us[2,l]){
    return(l_new[l])
  }
  return(l_current)
}

backward_pool <- function(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,c,delta,rw_const){
  # metropolis
  # sequential update
  X_pool <- array(0,dim=c(T,L,dim))
  
  # start
  es <- runif(L,0.1,0.4)
  es_2 <- es^2
  zs <- matrix(rnorm(L*dim),ncol=dim,nrow=L)
  Us <- runif(L)
  
  # sample the index for the current states 
  k <- sample(1:L,T,replace=TRUE)
  
  # place the first current state 
  X_pool[T,k[T],] <- X_current[T,]
  
  # transition down from k[T] to 1
  if(k[T]>1){
    this_X_current <- X_current[T,]
    for(l in (k[T]-1):1){
      # update latent states
      # X_new <- rw_update_x_oat_sys(this_X_current,X_pool,Y,t,l,l_current,rw_const,sigma_inv)
      # autoregressive update 
      X_new <- mu_init + sqrt(1-es_2[l])*(this_X_current-mu_init) + es[l]*sigma_init_L%*%zs[l,]
      # hastings ratio
      lnum <- sum(dpois(Y[T,],exp(c+delta*X_new),log=TRUE))
      ldenom <- sum(dpois(Y[T,],exp(c+delta*this_X_current),log=TRUE))
      hastings <- exp(lnum - ldenom)
      
      alpha <- min(1,hastings)
      if(alpha>Us[l]){
        this_X_current <- X_new
      } 
      # update l_current

      # set pool
      X_pool[T,l,] <- this_X_current
    } 
  }
  
  if(k[T]<L){
    this_X_current <- X_current[T,]
    for(l in (k[T]+1):L){
      # autoregressive update latent states
      X_new <- mu_init + sqrt(1-es_2[l])*(this_X_current-mu_init) + es[l]*sigma_init_L%*%zs[l,]
      # hastings ratio
      lnum <- sum(dpois(Y[T,],exp(c+delta*X_new),log=TRUE))
      ldenom <- sum(dpois(Y[T,],exp(c+delta*this_X_current),log=TRUE))
      hastings <- exp(lnum - ldenom)
      
      alpha <- min(1,hastings)
      if(alpha>Us[l]){
        this_X_current <- X_new
      } 
      
      # update l_current

      # set pool
      X_pool[T,l,] <- X_new
    }
  }
  
  acceptance_rate <- matrix(0,nrow=T-1,ncol=2)
  # then for t>1
  for(t in (T-1):1){
    #### stochastic initialisation of l 
    diff <- t(X_pool[t+1,,]) - as.vector(F %*% X_current[t,] )
    lprob_trans <- diag(-1/2*t(diff)%*%sigma_inv%*%diff)
    lprob_obs <- dpois(Y[t,],exp(c+delta*t(X_pool[t+1,,])),log=TRUE)
    lprob_obs <- apply(lprob_obs,MARGIN=2,FUN=sum)
    lprob <- lprob_obs + lprob_trans
    prob <- exp(lprob-max(lprob)); prob <- prob/sum(prob)
    l_original <- sample(1:L,1,prob=prob)
    ####

    this_X_pool <- matrix(logical(0),nrow=L,ncol=dim) # for this timestep only 
    this_X_pool[k[t],] <- X_current[t,]
    Us <- matrix(runif(L*2),nrow=2,ncol=L)
    l_new <- sample(1:L,L,replace=TRUE)

    # setup for normal update 
    es <- runif(L,0.1,0.4) # hardcoded
    zs <- matrix(rnorm(L*dim),nrow=L,ncol=dim) 
    last <- es*zs%*%sigma_U
    
    x_accept <- 0 # to measure acceptance rates 
    l_accept <- 0

    l_current <- l_original
    if(k[t]>1){
      this_X_current <- this_X_pool[k[t],]
      for(l in (k[t]-1):1){
        # update latent states
        X_new <- rw_update_x_oat_sys(this_X_current,X_pool,Y,t,l,l_current,rw_const,sigma_inv,Us)
        # X_new <- normal_update(this_X_current,X_pool,Y,t,l,l_current,rw_const,c,delta,sigma_inv,Us,last)

        if(all(this_X_current!=X_new)){
          x_accept <- x_accept + 1  
        }      
        this_X_current <- X_new

        # update l_current
        new_l <- uar_update_l(this_X_current,dim_index,X_pool,Y,t,l,l_current,c,delta,sigma_inv,l_new,Us)
        if(new_l!=l_current){
          l_accept <- l_accept+1
        }
        l_current <- new_l
        
        # set pool
        this_X_pool[l,] <- this_X_current
      } 
    }

    l_current <- l_original
    if(k[t]<L){
      this_X_current <- this_X_pool[k[t],]
      for(l in (k[t]+1):L){
        # this_X_current <- this_X_pool[l-1,]
        # update latent states
        X_new <- rw_update_x_oat_sys(this_X_current,X_pool,Y,t,l,l_current,rw_const,sigma_inv,Us)
        # X_new <- normal_update(this_X_current,X_pool,Y,t,l,l_current,rw_const,c,delta,sigma_inv,Us,last)

        if(all(this_X_current!=X_new)){
          x_accept <- x_accept + 1 
        }      
        this_X_current <- X_new
        
        # update l_current
        new_l <- uar_update_l(this_X_current,dim_index,X_pool,Y,t,l,l_current,c,delta,sigma_inv,l_new,Us)
        if(new_l!=l_current){
          l_accept <- l_accept+1
        }
        l_current <- new_l

        # set pool
        this_X_pool[l,] <- this_X_current
      } 
    }
    
    acceptance_rate[t,] <- c(x_accept,l_accept)/L
    X_pool[t,,] <- this_X_pool
    X_current[t,] <- this_X_current

  } 
  return(list(X_pool=X_pool,acceptance_rate=acceptance_rate))
}

forward_sampling <- function(X_pool,Y,L,T,dim,F,mu_init,sigma_init_inv,sigma_inv){
  # backward sampling
  X_new <- matrix(0,nrow=T,ncol=dim)
  # 2 parts: prior and observation model 
  # prior 
  diff_1 <- t(X_pool[1,,]) - mu_init
  lprob_prior <- diag(-1/2*t(diff_1)%*%sigma_init_inv%*%diff_1) # expect length L
  
  # emission 
  lprob_obs <- dpois(Y[1,],exp(c+delta*t(X_pool[1,,])),log=TRUE)
  lprob_obs <- apply(lprob_obs,MARGIN=2,FUN=sum)
  stopifnot(length(lprob_obs)==L)
  # total
  lprob <- lprob_obs + lprob_prior
  prob <- exp(lprob-max(lprob)); prob <- prob/sum(prob)
  # print(prob)
  X_new[1,] <- X_pool[1,sample(1:L,1,prob=prob),] 
  for(t in 2:T){
    diff <- t(X_pool[t,,]) - as.vector(F %*% X_new[t-1,] )
    lprob_trans <- diag(-1/2*t(diff)%*%sigma_inv%*%diff)
    lprob_obs <- dpois(Y[t,],exp(c+delta*t(X_pool[t,,])),log=TRUE)
    lprob_obs <- apply(lprob_obs,MARGIN=2,FUN=sum)
    lprob <- lprob_obs + lprob_trans

    prob <- exp(lprob-max(lprob)); prob <- prob/sum(prob)
    X_new[t,] <- X_pool[t,sample(1:L,1,prob=prob),]
  }
  return(X_new)
}

ehmmModel1_backward_RW <- function(ssm,N,L,init=NULL,seed=NULL){
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
  if((is.null(ssm$sigma_init_L)) | (is.null(ssm$sigma_init_U)) | (is.null(ssm$sigma_init_inv))){
    print('old ssm object, computing Choleskey and inverse for sigma_init')
    sigma_init_U <- chol(sigma_init)
    sigma_init_L <- t(sigma_init_U)
    sigma_init_inv <- chol2inv(sigma_init_U)
  } else {
    sigma_init_L <- ssm$sigma_init_L
    sigma_init_U <- ssm$sigma_init_U
    sigma_init_inv <- ssm$sigma_init_inv
  }
  
  if(is.null(ssm$sigma_inv)){
    print('old ssm object, computing inverse for sigma')
    sigma_inv <- chol2inv(sigma_U)
  } else {
    sigma_inv <- ssm$sigma_inv
  }
  
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
  rw_const <- runif(2*N)

  acceptance_rate <- array(0,dim=c(2*N,T-1,2)) # 2 for measuring only the autocorrelation and shift updates
  pb <- txtProgressBar(min=0,max=2*N,title="ehmm",style=3)

  # # debugging block for removing reversed sequence 
  # acceptance_rate <- array(0,dim=c(N,T-1,2))
  # pb <- txtProgressBar(min=0,max=N+1,title="pgbs",style=3)
  # X_sample <- array(0,dim=c(N+1,T,dim))
  
  for(i in seq(2,2*N,2)){
  # for(i in 2:(N+1)){
    # forward sequence
    pool_out <- backward_pool(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,c,delta,rw_const[i-1]) 

    X_pool <- pool_out$X_pool
    acceptance_rate[i-1,,] <- pool_out$acceptance_rate
    
    X_current <- forward_sampling(X_pool,Y,L,T,dim,F,mu_init,sigma_init_inv,sigma_inv)
    X_sample[i,,] <- X_current
    
    # reversed sequence
    pool_out <- backward_pool(X_current[seq(T,1,-1),],Y[seq(T,1,-1),],T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,c,delta,rw_const[i])

    X_pool <- pool_out$X_pool
    acceptance_rate[i,,] <- pool_out$acceptance_rate

    X_current[seq(T,1,-1),] <- forward_sampling(X_pool,Y[seq(T,1,-1),],L,T,dim,F,mu_init,sigma_init_inv,sigma_inv)
    X_sample[i+1,,] <- X_current
    
    setTxtProgressBar(pb, i)
  }
  return(list(X_sample=X_sample[-1,,],N=N,L=L,init=init,seed=seed,X_pool=X_pool,
              acceptance_rate=acceptance_rate))
}