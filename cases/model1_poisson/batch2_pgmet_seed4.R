setMKLthreads(2)
load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')
source('../../samplers/model1_pgmet.R')
N <- 70000
L <- 250
es <- c(0.2,0.8)
seed <- 4
filename <- paste('./data/batch2_pgbmet_seed',seed,'.RData',sep='')
print(filename)
system.time(
  pgmet_out <- pgmetModel1(ssm_poisson,N,L,es,seed=seed) 
)
save(pgmet_out,file=filename)