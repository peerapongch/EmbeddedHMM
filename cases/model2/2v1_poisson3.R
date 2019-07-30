rm(list=ls())
setMKLthreads(2)
load('./data/model2_poisson3_seed4.RData')
load('./data/model2_poisson3_seed4_env.RData')

source('../../samplers/model2_lambda_param_2.R')
N <- 100
L <- 80
N.mcmc.param <- 10
seed <- 1 
rw_scale <- c(0.1,0.1,0.1)
checkpoint.name <- './data/checkpoint_poisson3_2v1.RData'
system.time(
  lambda_out <- lambdaModel2_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,
  	seed=seed,rw_scale=rw_scale,checkpoint.name=checkpoint.name)
)

# save(lambda_out,file='./data/poisson3_lambda_param_2_v1.RData')

# load(file='./data/poisson3_lambda_param_2_v1.RData')
dim(lambda_out$param_sample)
library(coda)
x.param <- as.mcmc(lambda_out$param_sample)
plot(x.param)
plot(lambda_out$X_sample[,100,3])

# # load('./data/poisson3_lambda_param.RData')
# library(coda)
# x.param <- as.mcmc(lambda_out$param_sample)
# plot(x.param)
# 
# source('../../evaluation/function_mcmc_mu.R')
# param_mu <- apply(lambda_out$param_sample,MARGIN=2,FUN=mean)
# param_mu
# 
# plot(lambda_out$param_sample[sample(360001,3000),1],type='l')
# plot(lambda_out$X_sample[,100,3])
# # plot(ehmm_forward_out$X_sample[,100,3],ehmm_forward_out$X_sample[,100,2])
