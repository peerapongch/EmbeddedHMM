# evaluateAutocorrelationTime <- function(){}
load('./data/ssm_model1_poisson_1.R')

### different inits for the diagnostics ###
N <- 50000; e <- 0.5
init1 <- matrix(0,ncol=dim,nrow=T)
init2 <- matrix(1,ncol=dim,nrow=T)
init3 <- matrix(5,ncol=dim,nrow=T)
init4 <- matrix(10,ncol=dim,nrow=T)
init5 <- NULL

mcmcs <- list()
# mcmc <- mcmcGaussianSSM(10000,0.5,ssm_poisson,obs='Poisson')
# caution: ensure that the chain is run long enough so acf is almost zero for k>K
# profvis(mcmcGaussianSSM(N,e,ssm_poisson,obs='Poisson',init1))

### begin of the diagnostic function
mcmcs[[1]] <- mcmcGaussianSSM(N,e,ssm_poisson, obs='Poisson',init1)
mcmcs[[2]] <- mcmcGaussianSSM(N,e,ssm_poisson, obs='Poisson',init2)
mcmcs[[3]] <- mcmcGaussianSSM(N,e,ssm_poisson, obs='Poisson',init3)
mcmcs[[4]] <- mcmcGaussianSSM(N,e,ssm_poisson, obs='Poisson',init4)
mcmcs[[5]] <- mcmcGaussianSSM(N,e,ssm_poisson, obs='Poisson',init5)

save(mcmcs,file='./data/poisson_mcmcs_comparison_1.RData')
load(file='./data/poisson_mcmcs_comparison_2.RData')

# find the overall mean of each dimension at each time 
overall_mean <- matrix(0,nrow=T,ncol=dim)
remove <- round(N*0.1)
divisor <- N*length(mcmcs)
for(t in 1:T){
  for(j in 1:dim){
    sum <- 0
    for(i in 1:length(mcmcs)){
      sum <- sum + sum(mcmcs[[i]]$X_mcmc[,t,j])
    }
    overall_mean[t,j] <- sum/divisor
  }
}

# ### measure the actual sampling time
# # are all sampling times the same for different inits?
# # init1
# system.time(
#   replicate(5,mcmcGaussianSSM(N,e,ssm_poisson, obs='Poisson',init1))
# )/5 
# # output is time per 10000 samples
# # 0.0258 sec per sample 
# 
# # init 2
# t2 <- system.time(
#   mcmcGaussianSSM(N,e,ssm_poisson, obs='Poisson',init2)
# )
# # output is time per 10000 samples
# 
# # init 3
# t3 <- system.time(
#   mcmcGaussianSSM(N,e,ssm_poisson, obs='Poisson',init3)
# )
# # output is time per 10000 samples
# 
# # init 4
# t4 <- system.time(
#   mcmcGaussianSSM(N,e,ssm_poisson, obs='Poisson',init4)
# ) 
# # output is time per 10000 samples
# 
# # init 5
# t5 <- system.time(
#   mcmcGaussianSSM(N,e,ssm_poisson, obs='Poisson',init5)
# )
# # output is time per 10000 samples
# # all similar around 0.026

# find the autocorrelation time for each chain
ACTime <- function(mcmcs,T,dim,lag.max,tps){
  print(tps)
  # calculate overall mean 
  overall_mean <- matrix(0,nrow=T,ncol=dim)
  N <- mcmcs[[1]]$N
  remove <- round(N*0.1)
  N <- N-remove # reassign N  
  divisor <- N*length(mcmcs)
  for(t in 1:T){
    for(j in 1:dim){
      sum <- 0
      for(i in 1:length(mcmcs)){
        sum <- sum + sum(mcmcs[[i]]$X_mcmc[-(1:remove),t,j])
      }
      overall_mean[t,j] <- sum/divisor
    }
  }
  
  # calculate autocorrelation time
  ac_time <- array(0,dim=c(length(mcmcs),T,dim))
  pb <- txtProgressBar(min=0,max=length(mcmcs)*T*dim,style=3); prog <- 0
  for(i in 1:length(mcmcs)){
    X_mcmc <- mcmcs[[i]]$X_mcmc
    for(t in 1:T){
      for(j in 1:dim){
        setTxtProgressBar(pb, prog)
        prog <- prog + 1
        # calculate autocorrelation 
        chain <- X_mcmc[-(1:remove),t,j] - overall_mean[t,j]
        acfs <- acf(chain,demean=FALSE,lag.max=min(lag.max,length(chain)),type='correlation',plot=FALSE)
        rho_sum <- sum(acfs$acf)
        # adjust for time
        ac_time[i,t,j] <- (1+2*rho_sum)*tps
        print(ac_time[i,t,j])
      }
    }
  } 
  close(pb)
  return(ac_time)
  # return(list(ac_time=ac_time,overall_mean=overall_mean))
}

acf(mcmcs[[3]]$X_mcmc[,2,1],lag.max=50000)

actime_out <- ACTime(mcmcs,T,dim=1,lag.max=N-1,tps=0.026)
# save(actime_out,file='./data/actime_out_2.RData')
ac_time <- actime_out
# ac_time[1,1,1]
# plot(ac_time[1,,1],type='l')
# actime_out$overall_mean[1:10,1]
# overall_mean[1:10,1]


plot(ac_time[1,,1],type='l',ylim=c(0,40))
lines(ac_time[2,,1],col='red')
lines(ac_time[3,,1],col='blue')
lines(ac_time[4,,1],col='green')
lines(ac_time[5,,1],col='yellow')

mean_ac_time <- ac_time[1,,1]
for(i in 2:5){
  mean_ac_time <- mean_ac_time + ac_time[i,,1]
}
mean_ac_time <- mean_ac_time/5
plot(mean_ac_time,type='l')

autocov <- function(x,mean,k){
  n <- length(x)
  stopifnot(n>k)
  x_new <- x-mean
  sum <- 0 
  left <- x_new[1:(n-k)]
  right <- x_new[(k+1):n]
  gamma_k <- t(left) %*% right /n
  return(gamma_k[[1]])
}

X_mcmc <- mcmcs[[1]]$X_mcmc
x <- X_mcmc[,2,1]
o_mean <- overall_mean[2,1]
out <- rep(0,N-1)
for(i in 1:(N-1)){
  out[i] <- autocov(x,o_mean,i)
}
out_corr <- out/autocov(x,o_mean,0)
plot(out,type='l')
plot(out_corr,type='l')
tau <- 1+2*sum(out_corr)
tau_adj <- tau *tps

tau_adj <- rep(0,T)
for(t in 1:T){
  x <- X_mcmc[,t,1]
  o_mean <- overall_mean[t,1]
  out <- rep(0,N-1)
  for(i in 1:(N-1)){
    out[i] <- autocov(x,o_mean,i)
  }
  out_corr <- out/autocov(x,o_mean,0)
  # plot(out,type='l')
  # plot(out_corr,type='l')
  tau <- 1+2*sum(out_corr)
  tau_adj[t] <- tau *tps
}

acf(x,lag.max = 9999,type='covariance')

o_mean
mean(x)
# ac_time <- ACTime(mcmcs,T,dim,lag.max=5000)
# # find the autocorrelation for this set of chains 
# 
# tps <- 0.03 # time per sample 
# plot(ac_time[1,,1]*tps,type='l')

