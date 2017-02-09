t_meas_all <- c(seq(from=1/360, to=14/360, by=0.5/360), seq(from=15/360, to=30/360, by=1/360),
                seq(from=32/360, to=60/360, by=2/360), seq(from=67/360, to=360/360, by=7/360))


library(rstan)
fit1 <- readRDS("FITS/fit_1.rds")
f1 <- extract(fit1)
s1 <- summary(fit1)$summary
fit5 <- readRDS("FITS/fit_5.rds")
f5 <- extract(fit5)
s5 <- summary(fit5)$summary



Y1 <- colMeans(f1$CO2_flux_hat[,,1])
Y2 <- colMeans(f5$CO2_flux_hat[,,1])

plot(t_meas_all[seq(1,length(t_meas_all), 1)], Y1,type="l",ylim=c(0,0.04))
lines(t_meas_all[seq(1,length(t_meas_all), 5)], Y2, col="red")