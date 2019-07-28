setMKLthreads(2)
load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')

source('../../samplers/model1_lambda_param_RW_2.R')
N <- 9000
L <- 80
N.mcmc.param <- 20
seed <- 1 
rw_scale <- c(0.1,0.1,0.1,0.1)
checkpoint.name <- './data/checkpoint_poisson1_2v1.RData'
system.time(
  lambda_out <- lambdaModel1_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,
                                   seed=seed,
                                   rw_scale=rw_scale,
                                   checkpoint.name=checkpoint.name
  )
)

save(lambda_out,file='./data/poisson1_lambda_param_2_v1.RData')

library(coda)
x.param <- as.mcmc(lambda_out$param_sample)
plot(x.param)

plot(lambda_out$param_sample[,3])
