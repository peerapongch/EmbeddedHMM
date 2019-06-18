load('./data/ssm_model1_poisson_1.RData')
load('./data/ssm_model1_poisson_1_env.RData')
source('../../samplers/model1_pgbs.R')
library(Rfast)
library(profvis)

L <- 250
x_init <- matrix(0,nrow=T,ncol=dim) 
N <- 10
batch <- 4

system.time(
  pgbs_out <- pgbsModel1_batch(ssm_poisson,N,L,x_init,batch)
)
# save(pgbs_out,file='./data/pgbs_model1_poisson_1_L250_N2000.RData')

mci_mu <- matrix(0,nrow=T,ncol=dim)
for(t in 1:T){
  mci_mu[t,] <- mean(pgbs_out$X_sample[,t,])
}

plot(pgbs_out$X_sample[1,,1],col='grey',type='l',ylim=c(-20,6))
for(i in 1:N*batch){
  lines(pgbs_out$X_sample[i,,1],col='grey')
}
lines(ssm_poisson$X[,1],type='l',col='blue')
lines(mci_mu[,1],col='red')

source('../../evaluation/autocorrelation_time.R')
mcmcs <- list(pgbs_out)
ac_time <- ACTime(mcmcs,T,dim,1000,0.6)
plot(ac_time[1,,3],type='l')
