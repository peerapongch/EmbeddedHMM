# evaluateAutocorrelationTime <- function(){}
load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')
source('../../samplers/model1_mhmcmc.R')
# note: tps = 0.0364, lag.max=50 for 50,000 samples seems appropriate

### different inits for the diagnostics ###
N <- 10000; e <- 0.5
seeds <- seq(1,5,1)

mcmcs <- list()

### begin of the diagnostic function
for(i in 1:lengtH(seeds)){
  system.time(
    mcmcs[[i]] <- mcmcGaussianSSM(N,e,ssm_poisson, obs='Poisson',seeds[i])
  )
} # this takes 5 hours

# save(mcmcs,file='./data/poisson1_mcmc_comparison_set1.RData')
load(file='./data/poisson1_mcmcs_comparison_1.RData')

### re-form the set: 10 thinning ### 
source('../../evaluation/autocorrelation_time.R')
for(i in 1:length(mcmcs)){
  mcmcs[[i]]$X_sample <- mcmcs[[i]]$X_mcmc[seq(1,N,10),,]
}

# visualise the autocorrelation
acf(mcmcs[[3]]$X_sample[,2,1],lag.max=1000)

actime_out <- ACTime(mcmcs,T,dim=dim,lag.max=1000,tps=1)  

# save(actime_out,file='./data/actime_out_poisson1_comparison_set1_lag50.RData')
# load('./data/actime_out_poisson1_comparison_set1_lag50.RData')

ac_time <- actime_out$ac_time
plot(ac_time[1,,1],type='l',ylim=c(0,120))
lines(ac_time[2,,1],col='red')
lines(ac_time[3,,1],col='blue')
lines(ac_time[4,,1],col='green')
lines(ac_time[5,,1],col='yellow')

mean_ac_time <- ac_time[1,,1]
for(i in 2:5){
  mean_ac_time <- mean_ac_time + ac_time[i,,1]
}
mean_ac_time <- mean_ac_time/5
plot(mean_ac_time,type='l',col='red')
