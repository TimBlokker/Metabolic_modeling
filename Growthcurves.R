mgamrcm<-read.table("FECALIBACTERIUM PRAUSNITZII_mGAM_rcm", header = T)
time_temp<-read.table("time_temp", header=T)
as.character(time_temp[,1])
library(lubridate)
time<-as.numeric(lubridate::as.difftime(as.character(time_temp[,1])))
temp<-as.numeric(time_temp[,2])
rcm<-cbind(time,mgamrcm[,1:3])
mgam<-cbind(time,mgamrcm[,4:6])
head(rcm)
head(mgam)

library(ggplot2)
library(reshape2)
#theme_set(theme_bw)
ggplot2::theme_set(theme_classic())
rcm_high <- melt(rcm, id.vars="time")
head(rcm_high)
ggrcm<-ggplot(rcm_high, aes(x = time, y = value, col=variable)) + geom_point(alpha=0.7)
ggrcm

mgam_high <- melt(mgam, id.vars="time")
ggrcm<-ggplot(mgam_high, aes(x = time, y = value, col=variable)) + geom_point(alpha=0.7)
ggrcm


#install.packages("growthrates")
library(growthrates)
fit <- fit_easylinear(rcm$time, rcm$rcm2)
par(mfrow = c(1, 2))
plot(fit, log = "y")
plot(fit)


coef(fit)      # exponential growth parameters
mumax<-coef(fit)[3]
rsquared(fit)  # coefficient of determination (of log-transformed data)
