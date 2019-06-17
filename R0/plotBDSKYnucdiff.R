######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)


# get the names run log files
log <- list.files(path="./bdskynucdiff/out", pattern="*bdsky[_]subs[1234567890][.]log", full.names = TRUE)

for (i in seq(1, length(log))){
  # Read in the MASCOT *.logs
  t <- read.table(log[i], header=TRUE, sep="\t")
  # take a 10 % burn in
  t <- t[-seq(1,ceiling(length(t$Sample)/10)), ]
  # get all labels
  all_labels = labels(t)
  new.rate <- data.frame(s=rep(i, length(t$Sample)))
  new.sampling <- data.frame(s=rep(i, length(t$Sample)))
  
  # calculate the R0's in real space at all time points
  for (j in seq(1, length(all_labels[[2]]))){
    if (startsWith(all_labels[[2]][[j]],"absolute")){
      absoluteR0 = t[,all_labels[[2]][[j]]]
    }
    if (startsWith(all_labels[[2]][[j]],"log")){
      logeR0 = t[,all_labels[[2]][[j]]]
      new.rate[, all_labels[[2]][[j]]] <- exp(logeR0)*absoluteR0
    }
    if (startsWith(all_labels[[2]][[j]],"samplingProportion.1")){
      new.sampling[, "sampling.proportion"] <-  t[,all_labels[[2]][[j]]]
    }
  }
  
  if (i==1){
    rate=new.rate
    sampling = new.sampling
  }else{
    rate = rbind(rate,new.rate)
    sampling = rbind(sampling,new.sampling)
  }
  
}

# get all labels again
all_labels = labels(rate)
# calculate means, meadians upper and lower bounds
for (j in seq(2, length(all_labels[[2]]))){
  new.dat <- data.frame(t=j, mean=mean(rate[,all_labels[[2]][[j]]]), median=median(rate[,all_labels[[2]][[j]]]),
                       upper=quantile(rate[,all_labels[[2]][[j]]],0.975), lower=quantile(rate[,all_labels[[2]][[j]]],0.025))
  if (j==2){
    dat = new.dat
  }else{
    dat = rbind(dat, new.dat)
  }
}

# print the mean reproductive number to a file (for the simulations of persistance)
print_dat = dat
print_dat = print_dat[,"mean"]
print_dat$date = dat$date
write.table(rate, file="./../Persistence/r0.csv", sep=",")

# read in the weather data
weather.hour <- read.table("./../NonSequenceData/WeatherDataBasel.csv", header=TRUE, sep=";")
# and the school holidays
school <- read.table("./scoolholidays.csv", header=TRUE, sep=",")

# get the min max and mean values  over the day
for (i in seq(1, length(weather.hour$Year)/24)){
  init_point = (i-1)*24+1
  min.vals = weather.hour[init_point,]
  max.vals = weather.hour[init_point,]
  mean.vals = weather.hour[init_point,]
  for (j in seq(3, length(weather.hour))){
    min.vals[1,j] = min(weather.hour[init_point:(init_point+23), j])
    max.vals[1,j] = max(weather.hour[init_point:(init_point+23), j])
    mean.vals[1,j] = mean(weather.hour[init_point:(init_point+23), j])
  }
  min.vals$SchoolHolidays  = school[i, "School_Day"]
  max.vals$SchoolHolidays  = school[i, "School_Day"]
  mean.vals$SchoolHolidays  = school[i, "School_Day"]
  
  min.vals$bound = "min"
  max.vals$bound = "max"
  mean.vals$bound = "mean"
  
  min.vals$time = as.Date(paste(min.vals$Year,min.vals$Month,min.vals$Day,sep="-"))
  max.vals$time = as.Date(paste(min.vals$Year,min.vals$Month,min.vals$Day,sep="-"))
  mean.vals$time = as.Date(paste(min.vals$Year,min.vals$Month,min.vals$Day,sep="-"))
  
  new.vals = rbind(min.vals, max.vals, mean.vals)
  
  if (i==1){
    vals = new.vals
  }else{
    vals = rbind(vals, new.vals)
  }
}

dat$t = dat$t-max(dat$t)
dat$t = dat$t*2

# set the most recent sampled individual
mrsi = as.Date("2017-04-05")
dat$date = mrsi

for (i in seq(1, length(dat$t))){
  dat[i, "date"] = mrsi+dat[i,"t"]
}

# also read in case data
t_sub <- read.table("./../Subtypes/subtypes_distribution.csv", header = T, sep=",",check.names = F)

t_sub$date = as.Date(as.Date("2017-01-01"))
for (i in seq(1,length(t_sub$week))){
  t_sub[i,"date"] = as.Date(t_sub[i,"week"])
}

t_sub = t_sub[which(t_sub$clade=="A1"),]
case_date = t_sub[which(t_sub$all=="all"),]
case_date$numbers = case_date$numbers + t_sub[which(t_sub$all=="A1"),]$numbers


p_reff <- ggplot() +
  geom_bar(data=case_date, aes(x=date, y=numbers/50, fill=all), color=NA, stat="identity") 

p_reff <- p_reff + geom_line(data=dat, aes(x=date, y=mean, color="mean_Reff")) +
  geom_ribbon(data=dat, aes(x=date, ymin=upper, ymax=lower, fill="confidence interval"), alpha="0.2") +
  scale_y_continuous(limits=c(0, 2.6),sec.axis = sec_axis(~.*50, name = "weekly sequenced cases", breaks=c(0, 50,  100,  150,  200, 250))) +
  theme_minimal() + 
  xlab("") + 
  theme(axis.ticks.x=element_line(), legend.position = "none")+
  scale_x_date(date_labels = "%d %b %y", limits=c(as.Date("2016-11-22"), as.Date("2017-3-15")))+
  ylab("effective reproduction number") +
  scale_color_manual(name="",values = c("black"), labels=expression('mean R' [eff]))+
  scale_fill_manual(name="",values = c(rgb(0.8,0.8,0.8), "red"), labels=c("weekly cases",expression('95% interval')))
plot(p_reff)

weather_dat <- vals[which(vals$bound=="max"), c("time", "Temperature")]
weather_dat$TemperatureUpper = weather_dat$Temperature
weather_dat$TemperatureLower<- vals[which(vals$bound=="min"), "Temperature"]
weather_dat$TemperatureMean<- vals[which(vals$bound=="mean"), "Temperature"]

weather_dat$RelHumLower<- vals[which(vals$bound=="min"), "Relative_Humidity"]
weather_dat$RelHumUpper<- vals[which(vals$bound=="max"), "Relative_Humidity"]
weather_dat$RelHumMean<- vals[which(vals$bound=="mean"), "Relative_Humidity"]

weather_dat$SchoolHoliday<- vals[which(vals$bound=="mean"), "SchoolHolidays"]



p_wea <- ggplot() +
  geom_bar(data=case_date, aes(x=date, y=numbers/50, fill=all), color=NA, stat="identity") +
  geom_step(data = weather_dat, aes(x=time, y=SchoolHoliday*1.5+0.5, color="School Day")) +
  geom_line(data=dat, aes(x=date, y=mean, color="mean_Reff R_e")) +
  geom_ribbon(data=dat, aes(x=date, ymin=upper, ymax=lower, fill="confidence interval"), alpha="0.1") +
  geom_line(data = weather_dat, aes(x=time, y=RelHumMean/50, color="Humidity")) +
  geom_line(data = weather_dat, aes(x=time, y=TemperatureMean/20+1, color="Temperature")) +
  theme(axis.ticks.x=element_line())+
  scale_x_date(date_labels = "%d %b %y", limits=c(as.Date("2016-11-22"), as.Date("2017-3-15")))+
  scale_y_continuous(limits=c(0.5,2.3), breaks=c(0.5,1,1.5,2),labels=c(-10,0,10,20),
                     sec.axis = sec_axis(~.*50, name = "Relative humidity [%]", breaks=c(25,50, 75,100))) +
  theme_minimal() + 
  xlab("") + 
  ylab("Air temperature [°C]") +
  scale_color_manual(name="",values = c("mean_Reff R_e"="black", "Humidity"="#0072B2", "Temperature"="#D55E00","School Day" = "#009E73"), 
                     labels=c("Relative Humidity", expression('mean R' [eff] ), "Temperature",  "School Day"))+
    scale_fill_manual(name="",values = c(rgb(0.8,0.8,0.8), "red"), labels=c("weekly cases",expression('95% interval')))
plot(p_wea)

ggsave(plot=p_reff,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Global/bdsky_nucdiff.pdf",width=6, height=4)
ggsave(plot=p_wea,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Global/weather_nucdiff_legend.pdf",width=6, height=4)
p_wea = p_wea +   theme(legend.position = "none")
ggsave(plot=p_wea,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Global/weather_nucdiff.pdf",width=6, height=4)


## calculate cross-correlations between R_eff and the mean weather values
re.for.corr <- dat[which(dat$date>=as.Date("2016-12-01")),]
re.for.corr <- re.for.corr[which(re.for.corr$date<=as.Date("2017-2-01")),]
weather.for.corr <- vals[is.element(vals$time,re.for.corr$date),]
weather.for.corr <- weather.for.corr[which(weather.for.corr$bound=="mean"),]


weather.labels = labels(weather.for.corr)

library(Hmisc)
for (i in seq(6,length(weather.labels[[2]])-2)){
  weather.labels[[2]][i] 
  cor.vals = rcorr(re.for.corr$mean, weather.for.corr[,weather.labels[[2]][i]], type="pearson")
  new.corr.vals = data.frame(label=weather.labels[[2]][i], correlation=cor.vals$r["x","y"], p_value=cor.vals$P["x","y"])
  if (i==6){
    corr.vals = new.corr.vals
  }else{
    corr.vals = rbind(corr.vals, new.corr.vals)
  }
}
rownames(corr.vals)<-c()

write.table(corr.vals, "./Corr_Dez_to_Feb.tsv", sep="\t")


re.for.corr <- dat[which(dat$date>=as.Date("2016-11-01")),]
re.for.corr <- re.for.corr[which(re.for.corr$date<=as.Date("2017-3-01")),]
weather.for.corr <- vals[is.element(vals$time,re.for.corr$date),]
weather.for.corr <- weather.for.corr[which(weather.for.corr$bound=="mean"),]

weather.labels = labels(weather.for.corr)

library(Hmisc)
for (i in seq(6,length(weather.labels[[2]])-2)){
  weather.labels[[2]][i] 
  cor.vals = rcorr(re.for.corr$mean, weather.for.corr[,weather.labels[[2]][i]], type="pearson")
  new.corr.vals = data.frame(label=weather.labels[[2]][i], correlation=cor.vals$r["x","y"], p_value=cor.vals$P["x","y"])
  if (i==6){
    corr.vals2 = new.corr.vals
  }else{
    corr.vals2 = rbind(corr.vals2, new.corr.vals)
  }
}
rownames(corr.vals)<-c()

write.table(corr.vals, "./Corr_Nov_to_Mar.tsv", sep="\t")

corr.vals[,"p-value"] = corr.vals$p_value
corr.vals$p_value = c()
corr.vals2[,"p-value"] = corr.vals2$p_value
corr.vals2$p_value = c()

rownames(corr.vals) = corr.vals$label


library(grid)
library(gridExtra)
library(gtable)


title1 <- textGrob("2016-12-01 until 2017-2-01",gp=gpar(fontsize=10), just = "left")
title2 <- textGrob("2016-11-01 until 2017-3-01",gp=gpar(fontsize=10))

corr.vals$label = c()
corr.vals2$label = c()

tab1 = tableGrob(corr.vals)
tab2 = tableGrob(corr.vals2,rows = NULL)

padding <- unit(5,"mm")

table1 <- gtable_add_rows(
  tab1, 
  heights = grobHeight(title1) + padding,
  pos = 0)

table2 <- gtable_add_rows(
  tab2, 
  heights = grobHeight(title2) + padding,
  pos = 0)


table1 <- gtable_add_grob(
  table1, 
  title1, 
  1, 1, 1, ncol(table1))

table2 <- gtable_add_grob(
  table2, 
  title2, 
  1, 1, 1, ncol(table2))

p_tab = arrangeGrob(table1,table2,ncol=2)

ggsave(plot=p_tab,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Supplement/Correlation/WeatherCorr.all.pdf", width=13, height=6)




# plot the sampling proprotion
p.samp = ggplot(data = sampling, aes(x=sampling.proportion)) +
  geom_histogram(aes(y=..density..), binwidth=0.001, colour="black", fill="white") + xlab("proportion of sampled individuals")  +
  theme_minimal()+
  geom_vline(aes(xintercept=median(sampling.proportion, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)
plot(p.samp)
ggsave(plot=p.samp,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Clusters/sampling_proportion.pdf",width=3, height=3)
