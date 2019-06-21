load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')
source('../../samplers/model1_pgbs.R')

L <- 250
N <- 40000
seed <- 5

s <- system.time(
  pgbs_out <- pgbsModel1(ssm_poisson,N,L,seed=seed)
)
tpi <- s[[3]]/N
tps <- tpi/2

save(tpi,tps,pgbs_out,file=paste('./data/batch1_pgbs_seed',seed,'.RData',sep=''))