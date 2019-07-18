load('./data/model2_poisson3.RData')
load('./data/model2_poisson3_env.RData')

source('../../samplers/model2_ehmm_backward.R')
N <- 1000
L <- 80
seed <- 1
system.time(
  ehmm_backward_out <- ehmmModel2_backward(ssm_poisson,N,L,seed=seed)
)
# save(ehmm_out,file='./data/ehmm_backward_model1_poisson1_L50_N6000_AR_SH.RData')

plot(ehmm_backward_out$X_sample[,100,2])
