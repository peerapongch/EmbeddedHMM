rm(list=ls())
setMKLthreads(2)
load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')

source('../../samplers/model1_lambda_param_RW_2.R')
N <- 9000
L <- 80
N.mcmc.param <- 20
seed <- 1 
rw_scale <- c(0.01,0.01,0.01,0.01)
checkpoint.name <- './data/checkpoint_poisson1_2.RData'
filename <- './data/poisson1_lambda_param_2.RData'
print(filename)
system.time(
  lambda_out <- lambdaModel1_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,
                                   seed=seed,rw_scale=rw_scale,checkpoint.name=checkpoint.name
  )
)

save(lambda_out,file=filename)