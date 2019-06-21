load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')
source('../../samplers/model1_pgmet.R')

L <- 250
N <- 30000
seed <- 3

s <- system.time(
  pgmet_out <- pgmetModel1(ssm_poisson,N,L,seed=seed)
)
tpi <- s[[3]]/N
tps <- tpi/4 # 0.30

save(tpi,tps,pgmet_out,file=paste('./data/batch1_pgbs_seed',seed,'.RData',sep=''))
