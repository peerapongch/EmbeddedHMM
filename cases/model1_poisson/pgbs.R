load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')
source('../../samplers/model1_pgbs.R')
# source('../../evaluation/autocorrelation_time.R')

N <- 1000
L <- 250
seed <- 1

system.time(
  pgbs_out <- pgbsModel1(ssm_poisson,N,L,seed)
)

save(pgbs_out,file='./data/poisson1_pgbs_L250_N1000_seed1.RData')
# load(file='./data/pgbs_model1_poisson_1_L250_N2000.RData')

# pgbs_out <- list(X_sample=pgbs_out,N=2000)

mci_mu <- matrix(0,nrow=T,ncol=dim)
for(t in 1:T){
  mci_mu[t,] <- mean(pgbs_out$X_sample[,t,])
}

plot(pgbs_out$X_sample[1,,1],col='grey',type='l',ylim=c(-9,7))
for(i in 1:N){
  lines(pgbs_out$X_sample[i,,1],col='grey')
}
lines(ssm_poisson$X[,1],type='l',col='blue')
lines(mci_mu[,1],col='red')

mcmcs <- list(pgbs_out)
actime_out <- ACTime(mcmcs,T,dim,1000,0.26)
plot(actime_out$ac_time[1,,1],type='l')
