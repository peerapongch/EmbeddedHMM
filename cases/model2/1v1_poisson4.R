setMKLthreads(2)
load('./data/model2_poisson4.RData')
load('./data/model2_poisson4_env.RData')

source('../../samplers/model2_lambda_param.R')
N <- 9000
L <- 80
N.mcmc.param <- 20
seed <- 1 
rw_scale <- c(0.1,0.1,0.1)
checkpoint.name <- './data/checkpoint_poisson4_1v1.RData'
system.time(
  lambda_out <- lambdaModel2_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,seed=seed,rw_scale=rw_scale,checkpoint.name=checkpoint.name)
)

save(lambda_out,file='./data/poisson4_lambda_param_1_v1.RData')

# load('./data/poisson3_lambda_param_v2.RData')
# library(coda)
# x.param <- as.mcmc(lambda_out$param_sample)
# plot(x.param)

# plot(lambda_out$param_sample[,3])
# 
# plot(lambda_out$X_sample[,100,3])
# ssm_poisson$X[100,]

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
