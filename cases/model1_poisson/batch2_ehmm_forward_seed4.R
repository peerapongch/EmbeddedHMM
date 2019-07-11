setMKLthreads(2)
load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')
source('../../samplers/model1_ehmm_forward.R')
N <- 9000
L <- 50
seed <- 4
filename <- paste('./data/batch2_ehmm_forward_seed',seed,'.RData',sep='')
print(filename)
system.time(
  ehmm_forward_out <- ehmmModel1_forward(ssm_poisson,N,L,seed=seed) 
)
save(ehmm_forward_out,file=filename)