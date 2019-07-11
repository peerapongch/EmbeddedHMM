setMKLthreads(2)
load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')
source('../../samplers/model1_pgbs.R')
N <- 70000
L <- 250
seed <- 4
filename <- paste('./data/batch2_pgbs_seed',seed,'.RData',sep='')
print(filename)
system.time(
  pgbs_out <- pgbsModel1(ssm_poisson,N,L,seed=seed) 
)
save(pgbs_out,file=filename)