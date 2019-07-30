rm(list=ls())
setMKLthreads(2)
load('./data/model2_poisson3_seed4.RData')
load('./data/model2_poisson3_seed4_env.RData')

source('../../samplers/model2_lambda_param.R')
N <- 9000
L <- 80
N.mcmc.param <- 10
seed <- 1 
rw_scale <- c(0.01,0.01,0.01)
checkpoint.name <- './data/checkpoint_poisson3_1v1.RData'
system.time(
  lambda_out <- lambdaModel2_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,seed=seed,rw_scale=rw_scale,
                                   checkpoint.name=checkpoint.name)
)

# save(lambda_out,file='./data/poisson3_seed4_lambda_param_1_v1.RData')

load('./data/poisson3_seed4_lambda_param_1_v1.RData')
lambda_out$param_acceptance_rate
dim(lambda_out$param_sample)
x.param <- as.mcmc(lambda_out$param_sample[seq(1,180001,10),])
plot(x.param)
par(mfrow=c(3,1))
plot(lambda_out$param_sample[seq(1,180001,100),1],main='Version 1, scaling factor 0.01',ylab=bquote(phi))
plot(lambda_out$param_sample[seq(1,180001,100),2],main='Version 1, scaling factor 0.01',ylab=bquote(rho))
plot(lambda_out$param_sample[seq(1,180001,100),3],main='Version 1, scaling factor 0.01',ylab=bquote(sigma))

# dim(lambda_out$param_sample)
# colSums(lambda_out$param_sample[2000:18000,])/length(2000:18000)

plot(lambda_out$X_sample[,250,1])
ssm_poisson$X[250,]
# plot(lambda_out$param_sample[,3])
#
# plot(lambda_out$X_sample[,100,3])
# ssm_poisson$X[100,];l;

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
