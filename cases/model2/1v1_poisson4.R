setMKLthreads(2)
load('./data/model2_poisson4.RData')
load('./data/model2_poisson4_env.RData')

source('../../samplers/model2_lambda_param.R')
N <- 2000
L <- 80
N.mcmc.param <- 20
seed <- 1 
rw_scale <- c(0.01,0.01,0.01)
checkpoint.name <- './data/checkpoint_poisson4_1v1.RData'
system.time(
  lambda_out <- lambdaModel2_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,seed=seed,rw_scale=rw_scale,checkpoint.name=checkpoint.name)
)

save(lambda_out,file='./data/poisson4_lambda_param_1_v1_thatlongrunatmidnight.RData')

# load('./data/poisson3_lambda_param_v2.RData')
library(coda)
dim(lambda_out$param_sample)
x.param <- as.mcmc(lambda_out$param_sample[seq(1,80000,100),])
plot(x.param)

acf(lambda_out$param_sample[seq(1,80000,10),1],lag.max=2000)
acf(lambda_out$param_sample[seq(1,80000,100),1],lag.max=2000)

acf(lambda_out$param_sample[seq(1000,80000,10),1],lag.max=2000)
acf(lambda_out$param_sample[seq(2000,80000,10),1],lag.max=2000)

acf(lambda_out$param_sample[seq(1,80000,10),2],lag.max=2000)
acf(lambda_out$param_sample[seq(1,80000,100),2],lag.max=2000)

acf(lambda_out$param_sample[seq(1000,80000,10),2],lag.max=2000)
acf(lambda_out$param_sample[seq(2000,80000,10),2],lag.max=2000)

acf(lambda_out$param_sample[seq(1,80000,10),3],lag.max=2000)
acf(lambda_out$param_sample[seq(1,80000,100),3],lag.max=2000)

acf(lambda_out$param_sample[seq(1000,80000,10),3],lag.max=2000)
acf(lambda_out$param_sample[seq(2000,80000,10),3],lag.max=2000)
