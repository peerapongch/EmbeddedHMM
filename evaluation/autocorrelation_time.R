ACTime <- function(mcmcs,T,dim,lag.max,tps){
  # print(tps)
  # calculate overall mean 
  overall_mean <- matrix(0,nrow=T,ncol=dim)
  N <- dim(mcmcs[[1]]$X_sample)[1]
  remove <- round(N*0.1)
  N <- N-remove # reassign N  
  divisor <- N*length(mcmcs)
  for(t in 1:T){
    for(j in 1:dim){
      sums <- lapply(mcmcs,FUN=function(x){sum(x$X_sample[-(1:remove),t,j])})
      # sums <- 0 
      # for(i in 1:length(mcmcs)){
      #   sums <- sums + sum(mcmcs[[i]]$X_sample[-(1:remove),t,j])
      # }
      overall_mean[t,j] <- sum(as.numeric(as.character(sums)))/divisor
      # overall_mean[t,j] <- sum(sums)/divisor
      # print(overall_mean[t,j])
    }
  }
  
  # calculate autocorrelation time
  ac_time <- array(0,dim=c(length(mcmcs),T,dim))
  pb <- txtProgressBar(min=0,max=length(mcmcs)*T*dim,style=3); prog <- 0
  for(i in 1:length(mcmcs)){
    X_sample <- mcmcs[[i]]$X_sample
    for(t in 1:T){
      for(j in 1:dim){
        setTxtProgressBar(pb, prog)
        prog <- prog + 1
        # calculate autocorrelation 
        chain <- X_sample[-(1:remove),t,j] - overall_mean[t,j]
        # print(chain)
        acfs <- acf(chain,demean=FALSE,lag.max=min(lag.max,length(chain)),type='correlation',plot=FALSE)
        rho_sum <- sum(acfs$acf[-1])
        # adjust for time
        ac_time[i,t,j] <- (1+2*rho_sum)*tps
        # print(rho_sum)
      }
    }
  } 
  close(pb)
  # return(ac_time)
  return(list(ac_time=ac_time,overall_mean=overall_mean))
}
