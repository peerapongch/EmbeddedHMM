KFKS <- function(ssm){
  dim_x <- ssm$dim; T <- ssm$T; Y <- ssm$Y
  mu_init <- ssm$mu_init; sigma_init <- ssm$sigma_init
  F <- ssm$F; G <- ssm$G; Q <- ssm$Q; H <- ssm$H; R <- ssm$R
  
  sigma_predict <- array(0,dim=c(T,dim_x,dim_x))
  sigma_filtering <- array(0,dim=c(T,dim_x,dim_x))
  mu_predict <- matrix(0,ncol=dim_x,nrow=T)
  mu_filtering <- matrix(0,ncol=dim_x,nrow=T)
  v <- matrix(0,ncol=dim_x,nrow=T) # residuals
  K <- array(0,dim=c(T,dim_x,dim_x))
  S <- array(0,dim=c(T,dim_x,dim_x))
  
  # filtering
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
  
  # smoothing
  J <- array(0,dim=c((T-1),dim_x,dim_x)) # backward kalman gain 
  sigma_smoothing <- array(0,dim=c(T,dim_x,dim_x))
  mu_smoothing <- matrix(0,ncol=dim_x,nrow=T)
  
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
  
  return(list(mu_smoothing=mu_smoothing,sigma_smoothing=sigma_smoothing,J=J,
              mu_predict=mu_predict, sigma_predict=sigma_predict, mu_filtering=mu_filtering, sigma_filtering=sigma_filtering,
              v=v,K=K,S=S))
}
