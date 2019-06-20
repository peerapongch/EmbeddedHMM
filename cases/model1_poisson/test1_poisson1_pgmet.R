source('../../samplers/model1_pgmet.R')
load('./data/model1_poisson1.Rdata')
load('./data/model1_poisson1_env.Rdata')

L <- 250
N <- 20000
seed <- 1
es <- c(0.2,0.8)

system.time(
  pgmet_out <- pgmetModel1(ssm_poisson,N,L,es,seed=seed)
)/N/4

save(pgmet_out,file='./data/test1_poisson1_pgmet_L250_N20000_seed1.RData')