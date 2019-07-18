setMKLthreads(2)
load('./data/model2_poisson3.RData')
load('./data/model2_poisson3_env.RData')

source('../../samplers/model2_combined.R')
N <- 250
L <- 80
L_particles <- 5000
es <- c(0.3,1)
init <- matrix(1,nrow=T,ncol=dim)
N.mcmc <- 50
system.time(
  combined_out <- ehmmModel2_combined(ssm_poisson,N,L,L_particles,es,N.mcmc=N.mcmc,init=init)
)
save(combined_out,file='./data/poisson3_combined.RData')

# plot(combined_out$X_sample[,100,3])
# n <- 30
# plot(combined_out$X_sample[1:n,100,3],combined_out$X_sample[1:n,100,2],type='l')
# 
# combined_out$X_sample[1:n,100,3]
# combined_out$X_sample[1:n,100,2]
