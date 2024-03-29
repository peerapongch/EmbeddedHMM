autoregressive_update <- function(this_X_current,mu,Y,t,l,l_current,F,delta,sq_term,Us,last){
# autoregressive_update <- function(this_X_current,X_pool,Y,t,l,l_current,F,sigma_L,c,delta,es,es_2,zs,Us){
  # mu_current <- F %*% X_pool[t-1,l_current,]
  mu_current <- mu[l_current,]
  X_new <- mu_current + sq_term[l]*(this_X_current-mu_current)+last[l,]
  # hastings ratio
  lnum <- sum(dpois(Y[t,],delta*abs(X_new),log=TRUE))
  ldenom <- sum(dpois(Y[t,],delta*abs(this_X_current),log=TRUE))
  lhastings <- lnum - ldenom
  hastings <- exp(lhastings)
  alpha <- min(1,hastings)
  if(alpha>Us[1,l]){
    return(X_new)
  }
  return(this_X_current)
}

shift_update <- function(this_X_current,X_pool,Y,t,l,l_current,F,delta,Us,l_new){
  # shift update: update l then x
  # update of l UAR(1:L) done above since independent of l
  X_new <- this_X_current + F %*% (X_pool[t-1,l_new[l],]-X_pool[t-1,l_current,]) # faster abit
  lnum <- sum(dpois(Y[t,],delta*abs(X_new),log=TRUE))
  ldenom <- sum(dpois(Y[t,],delta*abs(this_X_current),log=TRUE))
  lhastings <- lnum - ldenom
  hastings <- exp(lhastings)
  alpha <- min(1,hastings)
  if(alpha>Us[2,l]){
    return(list(X_new=X_new,l_new=l_new[l]))
  } else {
    return(list(X_new=this_X_current,l_new=l_current))
  }
}

forward_pool <- function(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,delta){
  # metropolis
  # sequential update
  X_pool <- array(0,dim=c(T,L,dim))
  acceptance_rate <- matrix(0,nrow=T,ncol=2)
  ar_accept <- 0
  sh_accept <- 0
  
  # start
  es <- runif(L,0.1,0.4)
  es_2 <- es^2
  zs <- matrix(rnorm(L*dim),ncol=dim,nrow=L)
  Us <- runif(L)
  
  # sample the index for the current states 
  k <- sample(1:L,T,replace=TRUE)
  
  # place the first current state 
  X_pool[1,k[1],] <- X_current[1,]
  
  # reversed transition down from k[1] to 1
  if(k[1]>1){
    for(l in (k[1]-1):1){
      if(l %% 2 == 0){ # if even do usual
        # autoregressive update 
        X_new <- mu_init + sqrt(1-es_2[l])*(X_pool[1,l+1,]-mu_init) + es[l]*sigma_init_L%*%zs[l,]
        # hastings ratio
        lnum <- sum(dpois(Y[1,],delta*abs(X_new),log=TRUE))
        ldenom <- sum(dpois(Y[1,],delta*abs(X_pool[1,l+1,]),log=TRUE))
        lhastings <- lnum - ldenom
        hastings <- exp(lhastings)
        
        alpha <- min(1,hastings)
        if(alpha>Us[l]){
          X_pool[1,l,] <- X_new
          ar_accept <- ar_accept + 1
        } else {
          X_pool[1,l,] <- X_pool[1,l+1,]
        }
      } else { # if odd do flip 
        X_new <- -1*X_pool[1,l+1,]
        X_pool[1,l,] <- X_new
      }

    } 
  }
  
  # forward transition 
  if(k[1]<L){
    for(l in (k[1]+1):L){
      if(l %% 2 == 0){ # if even do flip 
        X_new <- -1*X_pool[1,l-1,]
        X_pool[1,l,] <- X_new
      } else { # if odd do usual
        # autoregressive update 
        X_new <- mu_init + sqrt(1-es_2[l])*(X_pool[1,l-1,]-mu_init) + es[l]*sigma_init_L%*%zs[l,]
        # hastings ratio
        lnum <- sum(dpois(Y[1,],delta*abs(X_new),log=TRUE))
        ldenom <- sum(dpois(Y[1,],delta*abs(X_pool[1,l-1,]),log=TRUE))
        lhastings <- lnum - ldenom
        hastings <- exp(lhastings)
        
        alpha <- min(1,hastings)
        if(alpha>Us[l]){
          X_pool[1,l,] <- X_new
          ar_accept <- ar_accept + 1
        } else {
          X_pool[1,l,] <- X_pool[1,l-1,]
        }
      }
    }
  }
  acceptance_rate[1,] <- c(ar_accept,sh_accept)/(L-1)
  
  # then for t>1
  for(t in 2:T){
    this_X_current <- X_current[t,]
    # stochastic initialisation of l 
    diff <- this_X_current - t(X_pool[t-1,,] %*% F) # for symmetric F
    l_lprob <- diag(-1/2*t(diff)%*%sigma_inv%*%diff)
    prob <- exp(l_lprob-max(l_lprob)); prob <- prob/sum(prob)
    l_original <- sample(1:L,1,prob=prob)
    
    es <- runif(L,0.1,0.4) # hardcoded
    zs <- matrix(rnorm(L*dim),nrow=L,ncol=dim) 
    Us <- matrix(runif(L*2),nrow=2,ncol=L)
    l_new <- sample(1:L,L,replace=TRUE)
    
    # precompute to save cost
    this_X_pool <- matrix(logical(0),nrow=L,ncol=dim) # for this timestep only 
    this_X_pool[k[t],] <- this_X_current
    mu <- X_pool[t-1,,] %*% F # symmetric F 
    sq_term <- sqrt(1-es^2)
    last <- es*zs%*%sigma_U

    ar_accept <- 0 # to measure acceptance rates 
    sh_accept <- 0 
    
    # reversed transition: shift then AR
    l_current <- l_original
    if(k[t]>1){
      this_X_current <- X_current[t,] # now only used to carry information between AR and shift update
      for(l in (k[t]-1):1){
        if(l %% 2 == 0){ # even then do usual
          # Shift then AR 
          # SHIFT HERE 
          shift_out <- shift_update(this_X_current,X_pool,Y,t,l,l_current,F,delta,Us,l_new)
          # count and then update
          if(all(this_X_current!=shift_out$X_new)){
            sh_accept <- sh_accept + 1
          }
          this_X_current <- shift_out$X_new
          l_current <- shift_out$l_new

          # AUTOREGRESSIVE HERE         
          X_new <- autoregressive_update(this_X_current,mu,Y,t,l,l_current,F,delta,sq_term,Us,last)
          # count and then update
          if(all(this_X_current!=X_new)){
            ar_accept <- ar_accept + 1 
          }
          this_X_current <- X_new 
        } else { # odd do flip
          this_X_current <- -1*this_X_current
          if(l_current %% 2 == 0){ # l is even, then -1
            l_current <- l_current - 1
          } else {
            l_current <- l_current + 1
          }
        }

        # set pool
        this_X_pool[l,] <- this_X_current
      } 
    }
    
    # forward transition: AR then shift 
    l_current <- l_original
    if(k[t]<L){
      this_X_current <- X_current[t,] # now only used to carry information between AR and shift update
      for(l in (k[t]+1):L){
        if(l %% 2 == 0){ # even do flip
          this_X_current <- -1*this_X_current
          if(l_current %% 2 == 0){ # l is even, then -1
            l_current <- l_current - 1
          } else {
            l_current <- l_current + 1
          }
        } else { # odd do usual
          # AR then shift
          # AUTOREGRESSIVE HERE 
          X_new <- autoregressive_update(this_X_current,mu,Y,t,l,l_current,F,delta,sq_term,Us,last)
          # count and update
          if(all(this_X_current!=X_new)){
            ar_accept <- ar_accept + 1 
          }
          this_X_current <- X_new 

          shift_out <- shift_update(this_X_current,X_pool,Y,t,l,l_current,F,delta,Us,l_new)
          # count and then update
          if(all(this_X_current!=shift_out$X_new)){
            sh_accept <- sh_accept + 1
          }
          this_X_current <- shift_out$X_new
          l_current <- shift_out$l_new
        }
        # set pool
        this_X_pool[l,] <- this_X_current
      } 
    }
    
    acceptance_rate[t,] <- c(ar_accept,sh_accept)/(L-1)
    X_pool[t,,] <- this_X_pool
    X_current[t,] <- this_X_current
  } 
  return(list(X_pool=X_pool,acceptance_rate=acceptance_rate))
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

makeSigma <- function(rho,dim){
  sigma <- matrix(rho,ncol=dim,nrow=dim)
  diag(sigma) <- 1
  return(sigma)
}

makeSigma_init <- function(rho,phi){
  v <- 1/sqrt((1-phi^2))
  sigma_init <- rho * v %*% t(v)
  diag(sigma_init) <- diag(sigma_init)/rho
  return(sigma_init)
}

param_rw_lprob <- function(param_formed,X_current,Y,T){
  # form
  delta <- param_formed$delta
  F <- param_formed$F
  phi <- diag(F)
  sigma_init <- param_formed$sigma_init
  sigma <- param_formed$sigma
  sigma_init_inv <- param_formed$sigma_init_inv
  sigma_inv <- param_formed$sigma_inv

  X_current_t <- t(X_current)
  mu <- phi * X_current_t
  diff <- X_current_t[,2:T] - mu[,1:(T-1)]
  # p(x_1:T|param)
  lprob1 <- -1/2*(log(det(sigma_init)) + (T-1)*log(det(sigma)) + t(X_current[1,])%*%sigma_init_inv%*%X_current[1,] + sum(diag(t(diff)%*%sigma_inv%*%diff)))
  # print("optimised")
  # print(lprob1)

  # # might this be wrong? --> try the unoptimised version
  # lprob1 <- dmvnorm(X_current[1,],c(0,0,0),sigma_init,log=TRUE)
  # print(lprob1)
  # for(t in 2:500){
  #   lprob1 <- lprob1 + dmvnorm(X_current[t,],mu[,t-1],sigma,log=TRUE)
  # }
  # print("unoptimised")
  # print(lprob1)

  # p(y_1:T|x_1:T,param)
  lambda <- as.vector(delta*abs(t(X_current)))
  lprob2 <- sum(dpois(as.vector(t(Y)),lambda,log=TRUE))

  return(lprob1+lprob2)
}

param_update_rw <- function(param_current,param_formed_current,param_lprob_current,X_current,Y,rw_scale,dim){
  # propose
  # param_new <- param_current + 2*rw_scale*(runif(3)-1/2)
  # param_new <- param_current + rw_scale*rnorm(3)
  
  param_new <- param_current + rw_scale*rnorm(3)
  param_new[3] <- c(0.8)
  check_domain_phi <- (param_new[1] > 0) && (param_new[1] < 1)
  check_domain_rho <- (param_new[2] > 0) && (param_new[2] < 1)
  # check_domain_delta <- (param_new[3] > 0)
  
  if(check_domain_phi && check_domain_rho){
    # form 
    F <- diag(rep(param_new[1],dim))
    delta_new <- rep(param_new[3],dim)
    sigma_init <- makeSigma_init(param_new[2],rep(param_new[1],dim))
    sigma <- makeSigma(param_new[2],dim)

    sigma_init_U <- chol(sigma_init)
    sigma_init_L <- t(sigma_init_U)
    sigma_init_inv <- chol2inv(sigma_init_U)
    sigma_U <- chol(sigma)
    sigma_L <- t(sigma_U)
    sigma_inv <- chol2inv(sigma_U)
    param_formed_new <- list(F=F,sigma_init=sigma_init,sigma=sigma,sigma_init_inv=sigma_init_inv,sigma_inv=sigma_inv,sigma_init_U = sigma_init_U, sigma_init_L = sigma_init_L, sigma_U = sigma_U, sigma_L = sigma_L,
      delta=delta_new)
    
    # compute prob
    param_lprob_new <- param_rw_lprob(param_formed_new,X_current,Y,T)
    # print('hastings old/new:')
    # print(param_lprob_current[1])
    # print(param_lprob_new[1])
    # print('---')
    if(log(runif(1))<(param_lprob_new-param_lprob_current)){
      # print('ACCEPTED!!')
      param_current <- param_new
      # print(param_current)
      param_lprob_current <- param_lprob_new
      param_formed_current <- param_formed_new
    }
  } # else automatic out and do nothing but return the current values

  return(list(param_new=param_current,param_lprob_new=param_lprob_current,param_formed_new=param_formed_current))
}

lambdaModel2_param <- function(ssm,N,L,N.mcmc.param=20,init=NULL,seed=NULL,rw_scale=c(0.5,0.5,0.5),rho.init=NULL,phi.init=NULL,delta.init=NULL,checkpoint.name=NULL){
  # note: one simplifying assumption is that we assume the knowledge of model 2 for parameterisation
  # sampling epsilon values instead
  # setup
  require(MASS)
  if(!is.null(seed)){
    set.seed(seed)
  }
  #poisson observation and gaussian latent process 
  Y <- ssm$Y
  dim <- ssm$dim
  T <- ssm$T
  mu_init <- ssm$mu_init # later include in the sampling step 

  if(is.null(phi.init)){
    phi.init <- runif(1)
  }
  if(is.null(rho.init)){
    rho.init <- runif(1)
  }
  if(is.null(delta.init)){
    delta.init <- 0.8
  }
  
  # parameters to be sampled
  param_sample <- matrix(logical(0),nrow=N.mcmc.param*2*N+1,ncol=3) # order is phi, rho; delta
  param_current <- c(phi.init,rho.init,delta.init)
  print(param_current)
  param_sample[1,] <- param_current

  # form
  delta <- delta.init
  F <- diag(rep(phi.init,dim))
  sigma_init <- makeSigma_init(rho.init,rep(phi.init,dim))
  sigma <- makeSigma(rho.init,dim)
  # then lots of precomputation 
  sigma_init_U <- chol(sigma_init)
  sigma_init_L <- t(sigma_init_U)
  sigma_init_inv <- chol2inv(sigma_init_U)
  sigma_U <- chol(sigma)
  sigma_L <- t(sigma_U)
  sigma_inv <- chol2inv(sigma_U)

  param_formed_current <- list(F=F,sigma_init=sigma_init,sigma=sigma,sigma_init_inv=sigma_init_inv,sigma_inv=sigma_inv,
    sigma_init_U=sigma_init_U, sigma_init_L=sigma_init_L, sigma_U=sigma_U, sigma_L=sigma_L, delta=delta)

  X_sample <- array(logical(0),dim=c(2*N+1,T,dim))
  if(is.null(init)){
    X_sample[1,,] <- mvrnorm(T,mu_init,sigma_init)
  } else {
    X_sample[1,,] <- init 
  }
  X_current <- X_sample[1,,]

  acceptance_rate <- array(0,dim=c(2*N,T,2)) # 2 for measuring only the autocorrelation and shift updates
  pb <- txtProgressBar(min=0,max=2*N,title="ehmm",style=3)
  param_acceptance_rate <- 0
  for(i in seq(2,2*N,2)){
    if(!is.null(checkpoint.name)){
      if(i == N){
        save(X_sample,param_sample,file=checkpoint.name)
      }
    }
    # do N.mcmc.param fix paramter updates between each latent sampling
    # latent: forward sequence
    pool_out <- forward_pool(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,delta)
    
    X_pool <- pool_out$X_pool
    acceptance_rate[i-1,,] <- pool_out$acceptance_rate
    
    X_current <- backward_sampling(X_pool,L,T,dim,F,sigma_inv)
    X_sample[i,,] <- X_current
    
    # parameter sampling 1
    param_lprob_current <- param_rw_lprob(param_formed_current,X_current,Y,T)
    start_index <- N.mcmc.param*(i-1) - (N.mcmc.param-1)
    for(j in 1:N.mcmc.param){
      param_update_out <- param_update_rw(param_current,param_formed_current,param_lprob_current,X_current,Y,rw_scale,dim)
      if(any(param_current != param_update_out$param_new)){
        param_acceptance_rate <- param_acceptance_rate + 1
      }
      param_formed_current <- param_update_out$param_formed_new
      param_current <- param_update_out$param_new
      param_lprob_current <- param_update_out$param_lprob_new
      param_sample[start_index+j,] <- param_current
    }
    
    # update parameters
    delta <- param_formed_current$delta
    F <- param_formed_current$F
    sigma_init <- param_formed_current$sigma_init
    sigma <- param_formed_current$sigma
    sigma_init_U <- param_formed_current$sigma_init_U
    sigma_init_L <- param_formed_current$sigma_init_L
    sigma_init_inv <- param_formed_current$sigma_init_inv
    sigma_U <- param_formed_current$sigma_U
    sigma_L <- param_formed_current$sigma_L
    sigma_inv <- param_formed_current$sigma_inv
    
    # latent: reversed sequence
    pool_out <- forward_pool(X_current[seq(T,1,-1),],Y[seq(T,1,-1),],T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,delta)

    X_pool <- pool_out$X_pool
    acceptance_rate[i,seq(T,1,-1),] <- pool_out$acceptance_rate

    X_current[seq(T,1,-1),] <- backward_sampling(X_pool,L,T,dim,F,sigma_inv)
    X_sample[i+1,,] <- X_current
    
    # parameter sampling 2 
    param_lprob_current <- param_rw_lprob(param_formed_current,X_current,Y,T)
    start_index <- N.mcmc.param*(i) - (N.mcmc.param-1)
    for(j in 1:N.mcmc.param){
      param_update_out <- param_update_rw(param_current,param_formed_current,param_lprob_current,X_current,Y,rw_scale,dim)
      if(any(param_current != param_update_out$param_new)){
        param_acceptance_rate <- param_acceptance_rate + 1
      }
      param_formed_current <- param_update_out$param_formed_new
      param_current <- param_update_out$param_new
      param_lprob_current <- param_update_out$param_lprob_new
      param_sample[start_index+j,] <- param_current
    }
    # update parameters
    delta <- param_formed_current$delta
    F <- param_formed_current$F
    sigma_init <- param_formed_current$sigma_init
    sigma <- param_formed_current$sigma
    sigma_init_U <- param_formed_current$sigma_init_U
    sigma_init_L <- param_formed_current$sigma_init_L
    sigma_init_inv <- param_formed_current$sigma_init_inv
    sigma_U <- param_formed_current$sigma_U
    sigma_L <- param_formed_current$sigma_L
    sigma_inv <- param_formed_current$sigma_inv
    
    setTxtProgressBar(pb, i)
  }
  param_acceptance_rate <- param_acceptance_rate/(N.mcmc.param*2*N)
  return(list(X_sample=X_sample[-1,,],param_sample=param_sample,N=N,L=L,init=init,
    seed=seed,X_pool=X_pool,acceptance_rate=acceptance_rate,param_acceptance_rate=param_acceptance_rate))
}
