setwd("C:/Users/peera/Documents/EmbeddedHMM/cases/model1_poisson/")
setMKLthreads(2)
load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')

source('../../samplers/model1_lambda_param_RW.R')
N <- 1000
L <- 80
N.mcmc.param <- 10
seed <- 1 
rw_scale <- c(0.01,0.01,0.01,0.01)
# rw_scale <- c(0.05,0.05,0.05,0.05)
checkpoint.name <- './data/checkpoint_poisson1_1v1.RData'
system.time(
  lambda_out <- lambdaModel1_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,
                                   seed=seed,
                                   rw_scale=rw_scale,
                                   checkpoint.name=checkpoint.name
                                   )
)

# save(lambda_out,file='./data/poisson1_lambda_param_1_v1_thatawesomemodel1run.RData')

library(coda)
lambda_out$param_acceptance_rate
dim(lambda_out$param_sample)
x.param <- as.mcmc(lambda_out$param_sample)
plot(x.param)

par(mfrow=c(2,1))
plot(lambda_out$param_sample[,1],main='scaling factor 0.01',ylab=bquote(phi))
plot(lambda_out$param_sample[,2],main='scaling factor 0.01',ylab=bquote(rho))
plot(lambda_out$param_sample[,3],main='scaling factor 0.01',ylab=bquote(sigma))
plot(lambda_out$param_sample[,3],main='scaling factor 0.01',ylab=bquote(c))

plot(lambda_out$param_sample[,3])

# all param, rerun unbounded c, 0.01
# good 
