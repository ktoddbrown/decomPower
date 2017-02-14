library(rstan)
# read all the fits
for (i in 1:6){
  for (j in 1:3) {
    assign(paste("f_",i,"_",j, sep=""),readRDS(paste("FITS/fit_",i,"_",j,".rds",sep=""))$fit)
    assign(paste("s_",i,"_",j, sep=""), 
    summary(readRDS(paste("FITS/fit_",i,"_",j,".rds",sep=""))$fit)$summary)
  }
}


params <- c("A2_g[1]")

X <- array(NA, c(18,2))
ind<-1
for (i in 1:6){
  for (j in 1:3) {
    X[ind,1] <- eval(as.name(paste("s_",i,"_",j, sep="")))[params,1]
    X[ind,2] <- eval(as.name(paste("s_",i,"_",j, sep="")))[params,3]
    ind <- ind+1
  }
}

coefplot(X[,1], X[,2],  CI=1)
points(c(1:18),X[,1]+0.1,col="red")