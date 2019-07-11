setMKLthreads(2)
load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')
source('../../samplers/model1_mhmcmc.R')
N <- 2000000
es <- c(0.2,0.8)
seed <- 5
filename <- paste('./data/batch2_mhmcmc_seed',seed,'.RData',sep='')
print(filename)
system.time(
  mhmcmc_out <- mhmcmcModel1(ssm_poisson,N,es,seed=seed) 
)
save(mhmcmc_out,file=filename)