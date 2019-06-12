# implements kalman smoother assuming the static parameters are known
# steps: 1. kalman filtering
# steps: find backward kalman gain to fint smoothing distribution

# init 
source('generate_model1_guassian.R')

# filtering: baby steps 
# filtering 1 
mu_1_0 <- F %*% mu_init
sigma_1_0 <- F %*% sigma_init %*% t(F) + G %*% Q %*% t(G)

S_1 <- H %*% sigma_1_0 %*% t(H) + R
K_1 <- sigma_1_0 %*% t(H) %*% solve(S_1)

v_1 <- Y[1,]-H %*% mu_1_0
mu_1_1 <- mu_1_0 + K_1 %*% v_1
sigma_1_1 <- (diag(1,dim_x) - K_1 %*% H) %*% sigma_1_0

##### filtering #####
# init: use array for implementation
sigma_predict <- array(0,dim=c(T,dim_x,dim_x))
sigma_filtering <- array(0,dim=c(T,dim_x,dim_x))
mu_predict <- matrix(0,ncol=dim_x,nrow=T)
mu_filtering <- matrix(0,ncol=dim_x,nrow=T)
v <- matrix(0,ncol=dim_x,nrow=T) # residuals
K <- array(0,dim=c(T,dim_x,dim_x))
S <- array(0,dim=c(T,dim_x,dim_x))

# begin: 
mu_predict[1,] <- F %*% mu_init
sigma_predict[1,,] <- F %*% sigma_init %*% t(F) + G %*% Q %*% t(G)
S[1,,] <- H %*% sigma_predict[1,,] %*% t(H) + R
K[1,,] <- sigma_predict[1,,] %*% t(H) %*% solve(S[1,,])
v[1,] <- Y[1,] - H %*% mu_predict[1,]
mu_filtering[1,] <- mu_predict[1,] + K[1,,] %*% v[1,]
sigma_filtering[1,,] <- (diag(1,dim_x) - K[1,,] %*% H) %*% sigma_predict[1,,]

for(t in 2:T){
  # predict
  mu_predict[t,] <- F %*% mu_filtering[t-1,]
  sigma_predict[t,,] <- F %*% sigma_filtering[t-1,,] %*% t(F) + G %*% Q %*% t(G)
  # update
  S[t,,] <- H %*% sigma_predict[t,,] %*% t(H) + R
  K[t,,] <- sigma_predict[t,,] %*% t(H) %*% solve(S[t,,])
  v[t,] <- Y[t,] - H %*% mu_predict[t,]
  mu_filtering[t,] <- mu_predict[t,] + K[t,,] %*% v[t,]
  sigma_filtering[t,,] <- (diag(1,dim_x) - K[t,,] %*% H) %*% sigma_predict[t,,] 
}

# filtering diagnostics 
plot(mu_filtering[,1],X[,2])
plot(X[,1],type='l')
lines(mu_filtering[,1],col='red')
points(Y[,1],col='blue')


# smoothing: baby steps
J_T1 <- sigma_filtering[T-1,,] %*% t(F) %*% solve(sigma_predict[T,,])
sigma_smooth_T1_T <- sigma_filtering[T-1,,] + J_T1 %*%(sigma_filtering[T,,]-sigma_predict[T,,])%*%t(J_T1)
mu_smooth_T1_T <- mu_filtering[T,] + J_T1 %*% (mu_filtering[T,] - mu_predict[T,]) # penultimate term: for this it is filtering, but for ther rest recurse on this

##### smoothing #####
# init: 
J <- array(0,dim=c((T-1),dim_x,dim_x)) # backward kalman gain 
sigma_smoothing <- array(0,dim=c(T,dim_x,dim_x))
mu_smoothing <- matrix(0,ncol=dim_x,nrow=T)

#begin: backward
mu_smoothing[T,] <- mu_filtering[T,]
sigma_smoothing[T,,] <- sigma_filtering[T,,]
for(t in (T-1):1){
  J[t,,] <- sigma_filtering[t,,] %*% t(F) %*% solve(sigma_predict[t+1,,])
  if(t==(T-1)){
    sigma_smoothing[t,,] <- sigma_filtering[t,,] + J[t,,] %*%(sigma_filtering[t+1,,]-sigma_predict[t+1,,])%*%t(J[t,,])
    mu_smoothing[t,] <- mu_filtering[t,] + J[t,,] %*% (mu_filtering[t+1,] - mu_predict[t+1,])
  } else {
    sigma_smoothing[t,,] <- sigma_filtering[t,,] + J[t,,] %*%(sigma_smoothing[t+1,,]-sigma_predict[t+1,,])%*%t(J[t,,])
    mu_smoothing[t,] <- mu_filtering[t,] + J[t,,] %*% (mu_smoothing[t+1,] - mu_predict[t+1,])
  }
}

# filtering + smoothing diagnostics 
plot(X[,1],type='l')
lines(mu_smoothing[,1],col='blue')
lines(mu_filtering[,1],col='red')


# complete: now package this and re implement the gaussian model so that the matrix values are encapsulated
