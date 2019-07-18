load('./data/model2_poisson3.RData')
load('./data/model2_poisson3_env.RData')

source('../../samplers/model2_combined.R')
N <- 1000
L <- 80
L_particles <- 500
es <- c(0.2,0.8)
seed <- 1
system.time(
  combined_out <- ehmmModel2_combined(ssm_poisson,N,L,L_particles,es,seed=seed)
)
# save(ehmm_out,file='./data/ehmm_backward_model1_poisson1_L50_N6000_AR_SH.RData')

plot(combined_out$X_sample[,100,3])
n <- 30
plot(combined_out$X_sample[1:n,100,3],combined_out$X_sample[1:n,100,2],type='l')

combined_out$X_sample[1:n,100,3]
combined_out$X_sample[1:n,100,2]
