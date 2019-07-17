setMKLthreads(2)
load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')
source('../../samplers/model1_ehmm_backward.R')
N <- 20000
L <- 50
elim <- c(0,0.1)
seed <- 1
filename <- paste('./data/batch2_ehmm_backward_seed',seed,'.RData',sep='')
print(filename)
system.time(
  ehmm_backward_out <- ehmmModel1_backward(ssm_poisson,N,L,seed=seed,elim=elim) 
)
save(ehmm_backward_out,file=filename)