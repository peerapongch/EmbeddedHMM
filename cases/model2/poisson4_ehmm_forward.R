load('./data/model2_poisson4.RData')
load('./data/model2_poisson4_env.RData')

source('../../samplers/model2_ehmm_forward.R')
N <- 9000
L <- 80
init <- matrix(1,nrow=T,ncol=dim)
system.time(
  ehmm_forward_out <- ehmmModel2_forward(ssm_poisson,N,L,init=init)
)
save(ehmm_forward_out,file='./data/poisson4_ehmm_forward.RData')

# plot(ehmm_forward_out$X_sample[,100,3],ehmm_forward_out$X_sample[,100,2])
