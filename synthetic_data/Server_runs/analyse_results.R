library(rstan)
library(arm)
# read all the fits
for (i in 1:7){
  for (j in 1) {
    assign(paste("f_",i,"_",j, sep=""),readRDS(paste("FITS/fit_",i,"_",j,".rds",sep=""))$fit)
    assign(paste("s_",i,"_",j, sep=""), 
    summary(readRDS(paste("FITS/fit_",i,"_",j,".rds",sep=""))$fit)$summary)
  }
}


params <- c("A1_g[2]")
X <- array(NA, c(7,2))
ind<-1
for (i in 1:7){
  for (j in 1) {
    X[ind,1] <- eval(as.name(paste("s_",i,"_",j, sep="")))[params,1]
    X[ind,2] <- eval(as.name(paste("s_",i,"_",j, sep="")))[params,3]
    ind <- ind+1
  }
}
downsample_set <- rep(c(1,2,4,5,10,20,50), each=1)
set_size <- 100/downsample_set
coefs <- data.frame(vars = sprintf("%02d",set_size))
coefs$est = X[,1];
coefs$se = X[,2];
ggplot(coefs, aes(vars, est)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=est - 1.96*se, ymax=est + 1.96*se), 
                lwd=1, colour="red", width=0) +
  geom_errorbar(aes(ymin=est - se, ymax=est + se), 
                lwd=2.5, colour="blue", width=0) +
  geom_point(size=4, pch=21, fill="yellow") +
  theme_bw()






fff<-readRDS("FITS/fit_4_1.rds")

library(rstan)
library(ggplot2);
fit <- extract(fff$fit);
num_rep<-5;
CO2_flux_hat <- array(0,c(fff$data.N_t, num_rep));
CO2_flux_lci <- array(NA,c(fff$data.N_t, num_rep));
CO2_flux_uci <- array(NA,c(fff$data.N_t, num_rep));
for (nr in 1:num_rep) {
  for (t in 1:fff$data.N_t) {
    CO2_flux_hat[t,nr] <- median(fit$CO2_flux_hat[,t,nr]);
    CO2_flux_lci[t,nr] <- quantile(fit$CO2_flux_hat[,t,nr], 0.025);
    CO2_flux_uci[t,nr] <-  quantile(fit$CO2_flux_hat[,t,nr], 0.975);
  }
}
df_plot <- data.frame(list(t_meas = fff$data.t_meas, 
                           CO2_flux_meas = fff$data.CO2_flux,
                           CO2_flux_hat = CO2_flux_hat, CO2_flux_lci = CO2_flux_lci, 
                           CO2_flux_uci = CO2_flux_uci));

ggplot(df_plot, aes(x = 360*t_meas)) +
  geom_ribbon(aes(ymin = CO2_flux_lci[,1],
                  ymax = CO2_flux_uci[,1]),
              fill="lightyellow") +
  geom_line(aes(y=CO2_flux_hat[,1]),colour="darkred") +
  geom_ribbon(aes(ymin = CO2_flux_lci[,2],
                  ymax = CO2_flux_uci[,2]),
              fill="lightyellow") +
  geom_line(aes(y=CO2_flux_hat[,2]),colour="darkred") +
  geom_ribbon(aes(ymin = CO2_flux_lci[,3],
                  ymax = CO2_flux_uci[,3]),
              fill="lightyellow") +
  geom_line(aes(y=CO2_flux_hat[,3]),colour="darkred") +
  geom_ribbon(aes(ymin = CO2_flux_lci[,4],
                  ymax = CO2_flux_uci[,4]),
              fill="lightyellow") +
  geom_line(aes(y=CO2_flux_hat[,4]),colour="darkred") +
  geom_ribbon(aes(ymin = CO2_flux_lci[,5],
                  ymax = CO2_flux_uci[,5]),
              fill="lightyellow") +
  geom_line(aes(y=CO2_flux_hat[,5]),colour="darkred") +
  geom_line(aes(y=fff$data.CO2_flux[,1]),colour="darkblue") +
  geom_line(aes(y=fff$data.CO2_flux[,2]),colour="darkblue") +
  geom_line(aes(y=fff$data.CO2_flux[,3]),colour="darkblue") +
  geom_line(aes(y=fff$data.CO2_flux[,4]),colour="darkblue") +
  geom_line(aes(y=fff$data.CO2_flux[,5]),colour="darkblue") +
  labs(x="Time (days)", 
       y="CO2 flux") +
  ggtitle("Soil Incubation: Century Model (hierarchical-nonstiff),  \n data (blue), estimate median (red),  95% interval (yellow)");
