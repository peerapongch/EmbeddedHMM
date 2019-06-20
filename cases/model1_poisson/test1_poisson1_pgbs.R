load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')
source('../../samplers/model1_pgbs.R')
source('../../evaluation/autocorrelation_time.R')

L <- 250
N <- 20000
seed <- 1

system.time(
  pgbs_out <- pgbsModel1(ssm_poisson,N,L,seed)
)

save(pgbs_out,file='./data/test1_poisson1_pgbs_L250_N20000_seed1.RData')
