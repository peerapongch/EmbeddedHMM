rm(list=ls())
setMKLthreads(2)
load('./data/model2_poisson3.RData')
load('./data/model2_poisson3_env.RData')

source('../../samplers/model2_lambda_param_3.R')
N <- 9000
L <- 80
N.mcmc.param <- 20
seed <- 1 
rw_scale <- c(0.01,0.01,0.01)
checkpoint.name <- './data/checkpoint_poisson3_3.RData'
filename <- './data/poisson3_lambda_param_3.RData'
print(filename)
system.time(
  lambda_out <- lambdaModel2_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,
                                   seed=seed,rw_scale=rw_scale,checkpoint.name=checkpoint.name
  )
)

save(lambda_out,file=filename)