install.packages('KFKSDS')
install.packages('stsm')
library(KFKSDS)
library(stsm)

# local level plus seasonal model with arbitrary parameter values
# for the 'JohnsonJohnson' time series
m <- stsm::stsm.model(model = "llm+seas", y = JohnsonJohnson,
                      pars = c("var1" = 2, "var2" = 15, "var3" = 30))
ss <- stsm::char2numeric(m)

# run the Kalman filter
kf <- KF(m@y, ss)
plot(kf$a.upd[,1:2], main = "filtered state vector")
# 'KF.C' is a faster version that returns only the
# value of the negative of the likelihood function
kfc <- KF.C(m@y, ss)
all.equal(kf$mll, kfc)

# compute also derivative terms used below
kfd <- KF.deriv(m@y, ss)
all.equal(kfc, kfd$mll)
kfdc <- KF.deriv.C(m@y, ss, return.all = TRUE)
all.equal(kf$mll, kfdc$mll)


# then smoother 
ks <- KS(m@y, ss, kf)
ks
plot(ks$ahat[,1:2], main = "smoothed state vector")

kfd <- KF.deriv(m@y, ss)
ksd <- KS.deriv(m@y, ss, kfd)
all.equal(ks$ahat, ksd$ahat)

