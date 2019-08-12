######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(grid)
library(gridExtra)

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


# list with the plots
p.samp = list()
p_reff = list()
p_wea = list()
p_tab = list()


for (sub in seq(1,10)){
  print(sub)
  # get the names run log files
  log <- list.files(path="./bdskynucdiff/out", pattern=paste("*bdsky[_]subs", sub ,"[.]log", sep=""), full.names = TRUE)
  
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
  print_dat = data.frame(mean=print_dat[,"mean"])
  print_dat$date = dat$date
  write.table(rate, file="./../Persistence/r0.csv", sep=",")
  
  # read in the weather data
  weather.hour <- read.table("./../NonSequenceData/WeatherDataBasel.csv", header=TRUE, sep=";")
  # and the school holidays
  school <- read.table("./scoolholidays.csv", header=TRUE, sep=",")
  
  # get the min max data over the day
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
  # Read in the percentiles
  t_sub <- read.table("./../Subtypes/subtypes_distribution.csv", header = T, sep=",",check.names = F)
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # plot the percentiles as heatmaps
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  t_sub$date = as.Date(as.Date("2017-01-01"))
  for (i in seq(1,length(t_sub$week))){
    t_sub[i,"date"] = as.Date(t_sub[i,"week"])
  }
  
  t_sub = t_sub[which(t_sub$clade=="A1"),]
  case_date = t_sub[which(t_sub$all=="all"),]
  case_date$numbers = case_date$numbers + t_sub[which(t_sub$all=="A1"),]$numbers
  
  
  p_reff[[sub]] <- ggplot() +
    geom_bar(data=case_date, aes(x=date, y=numbers/50, fill=all), color=NA, stat="identity") 
  
  p_reff[[sub]] <- p_reff[[sub]] + geom_line(data=dat, aes(x=date, y=mean, color="mean_Reff")) +
    geom_ribbon(data=dat, aes(x=date, ymin=upper, ymax=lower, fill="confidence interval"), alpha="0.2") +
    scale_y_continuous(limits=c(0, 2.6),sec.axis = sec_axis(~.*50, name = "weekly sequenced cases", breaks=c(0,50,  100,  150,  200, 250))) +
    theme_minimal() + 
    xlab("") + 
    theme(axis.ticks.x=element_line())+
    scale_x_date(date_labels = "%d %b %y", limits=c(as.Date("2016-11-22"), as.Date("2017-3-15")))+
    ylab("effective reproduction number") +
    scale_color_manual(name="",values = c("black"), labels=expression('mean R' [eff]))+
    scale_fill_manual(name="",values = c(rgb(0.8,0.8,0.8), "red"), labels=c("weekly cases",expression('95% interval'))) +
    theme(legend.position = "none")

  weather_dat <- vals[which(vals$bound=="max"), c("time", "Temperature")]
  weather_dat$TemperatureUpper = weather_dat$Temperature
  weather_dat$TemperatureLower<- vals[which(vals$bound=="min"), "Temperature"]
  weather_dat$TemperatureMean<- vals[which(vals$bound=="mean"), "Temperature"]
  
  weather_dat$RelHumLower<- vals[which(vals$bound=="min"), "Relative_Humidity"]
  weather_dat$RelHumUpper<- vals[which(vals$bound=="max"), "Relative_Humidity"]
  weather_dat$RelHumMean<- vals[which(vals$bound=="mean"), "Relative_Humidity"]
  
  weather_dat$SchoolHoliday<- vals[which(vals$bound=="mean"), "SchoolHolidays"]
  
  p_wea[[sub]] <- ggplot() +
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
    scale_color_manual(name="",values = c("mean_Reff R_e"="black", "Humidity"="#0072B2","School Day" = "#009E73", "Temperature"="#D55E00"), 
                       labels=c("Relative Humidity", expression('mean R' [eff] ),  "School Day", "Temperature"))+
    scale_fill_manual(name="",values = c(rgb(0.8,0.8,0.8), "red"), labels=c("weekly cases",expression('95% interval')))+
    theme(legend.position = "none")
  

  startDate = c(as.Date("2016-11-01"), as.Date("2016-12-01"))
  endDate = c(as.Date("2017-3-01"), as.Date("2017-2-01"))
  for (c in seq(1,length(startDate))){
    print(c)
    
    ## calculate correlations between R_eff and the mean weather values
    re.for.corr <- dat[which(dat$date>=startDate[[c]]),]
    re.for.corr <- re.for.corr[which(re.for.corr$date<=endDate[[c]]),]
    
    
    weather.for.corr <- vals[is.element(vals$time,re.for.corr$date),]
    weather.for.corr <- weather.for.corr[which(weather.for.corr$bound=="mean"),]
    
    # Average values over 2 days ---------------------------------
    
    
    # average the weather over 2 and 4 days
    remove(weather.for.corr.2days)
    for (i in seq(1,length(re.for.corr$dat), 1) ){
      tmp.vals1 = vals[is.element(vals$time,re.for.corr[i,"date"]),]
      tmp.vals1 = tmp.vals1[which(tmp.vals1$bound=="mean"),]
      
      tmp.vals2 = vals[is.element(vals$time,re.for.corr[i,"date"]+1),]
      tmp.vals2 = tmp.vals2[which(tmp.vals2$bound=="mean"),]
      
      # take the average
      new.data = tmp.vals1
      for (j in seq(1,length(new.data)-2)){
        new.data[1,j] = (new.data[1,j] +tmp.vals2[1,j])/2
      }
      if (i==1){
        weather.for.corr.2days = new.data
      }else{
        weather.for.corr.2days = rbind(weather.for.corr.2days,new.data)
      }
    }
    
    # Average values over 4 days ---------------------------------
    
    remove(weather.for.corr.4days)
    remove(re.for.corr.4days)
    
    for (i in seq(1,length(re.for.corr$dat), 2) ){
      tmp.vals1 = vals[is.element(vals$time,re.for.corr[i,"date"]),]
      tmp.vals1 = tmp.vals1[which(tmp.vals1$bound=="mean"),]
      
      tmp.vals2 = vals[is.element(vals$time,re.for.corr[i,"date"]+1),]
      tmp.vals2 = tmp.vals2[which(tmp.vals2$bound=="mean"),]
      
      tmp.vals3 = vals[is.element(vals$time,re.for.corr[i,"date"]+2),]
      tmp.vals3 = tmp.vals3[which(tmp.vals3$bound=="mean"),]
      
      tmp.vals4 = vals[is.element(vals$time,re.for.corr[i,"date"]+3),]
      tmp.vals4 = tmp.vals4[which(tmp.vals4$bound=="mean"),]
      
      
      # take the average
      new.data = tmp.vals1
      new.re.for.corr.4days = re.for.corr[i,]
      for (j in seq(1,length(new.data)-2)){
        new.data[1,j] = (new.data[1,j] + tmp.vals2[1,j] + tmp.vals3[1,j] + tmp.vals4[1,j])/4
      }
      for (j in seq(1,length(re.for.corr)-1)){
        new.re.for.corr.4days[1,j] = (new.re.for.corr.4days[1,j]+ re.for.corr[i+1,j])/2
      }
      if (i==1){
        weather.for.corr.4days = new.data
        re.for.corr.4days = new.re.for.corr.4days
      }else{
        weather.for.corr.4days = rbind(weather.for.corr.4days,new.data)
        re.for.corr.4days = rbind(re.for.corr.4days,new.re.for.corr.4days)
      }
    }
    
    # Average values over 6 days ---------------------------------
    
    
    remove(weather.for.corr.6days)
    remove(re.for.corr.6days)
    
    for (i in seq(1,length(re.for.corr$dat), 3) ){
      tmp.vals1 = vals[is.element(vals$time,re.for.corr[i,"date"]),]
      tmp.vals1 = tmp.vals1[which(tmp.vals1$bound=="mean"),]
      
      tmp.vals2 = vals[is.element(vals$time,re.for.corr[i,"date"]+1),]
      tmp.vals2 = tmp.vals2[which(tmp.vals2$bound=="mean"),]
      
      tmp.vals3 = vals[is.element(vals$time,re.for.corr[i,"date"]+2),]
      tmp.vals3 = tmp.vals3[which(tmp.vals3$bound=="mean"),]
      
      tmp.vals4 = vals[is.element(vals$time,re.for.corr[i,"date"]+3),]
      tmp.vals4 = tmp.vals4[which(tmp.vals4$bound=="mean"),]
      
      tmp.vals5 = vals[is.element(vals$time,re.for.corr[i,"date"]+4),]
      tmp.vals5 = tmp.vals5[which(tmp.vals5$bound=="mean"),]
      
      tmp.vals6 = vals[is.element(vals$time,re.for.corr[i,"date"]+5),]
      tmp.vals6 = tmp.vals6[which(tmp.vals6$bound=="mean"),]
      
      
      # take the average
      new.data = tmp.vals1
      new.re.for.corr.6days = re.for.corr[i,]
      for (j in seq(1,length(new.data)-2)){
        new.data[1,j] = (new.data[1,j] + tmp.vals2[1,j] + tmp.vals3[1,j] + tmp.vals4[1,j]  + tmp.vals5[1,j]  + tmp.vals6[1,j])/6
      }
      for (j in seq(1,length(re.for.corr)-1)){
        new.re.for.corr.6days[1,j] = (new.re.for.corr.6days[1,j] + re.for.corr[i+1,j] + re.for.corr[i+2,j])/3
      }
      if (i==1){
        weather.for.corr.6days = new.data
        re.for.corr.6days = new.re.for.corr.6days
      }else{
        weather.for.corr.6days = rbind(weather.for.corr.6days,new.data)
        re.for.corr.6days = rbind(re.for.corr.6days,new.re.for.corr.6days)
      }
    }
    
    # Average values over 8 days ---------------------------------
    
    
    remove(weather.for.corr.8days)
    remove(re.for.corr.8days)
    
    for (i in seq(1,length(re.for.corr$dat), 4) ){
      tmp.vals1 = vals[is.element(vals$time,re.for.corr[i,"date"]),]
      tmp.vals1 = tmp.vals1[which(tmp.vals1$bound=="mean"),]
      
      tmp.vals2 = vals[is.element(vals$time,re.for.corr[i,"date"]+1),]
      tmp.vals2 = tmp.vals2[which(tmp.vals2$bound=="mean"),]
      
      tmp.vals3 = vals[is.element(vals$time,re.for.corr[i,"date"]+2),]
      tmp.vals3 = tmp.vals3[which(tmp.vals3$bound=="mean"),]
      
      tmp.vals4 = vals[is.element(vals$time,re.for.corr[i,"date"]+3),]
      tmp.vals4 = tmp.vals4[which(tmp.vals4$bound=="mean"),]
      
      tmp.vals5 = vals[is.element(vals$time,re.for.corr[i,"date"]+4),]
      tmp.vals5 = tmp.vals5[which(tmp.vals5$bound=="mean"),]
      
      tmp.vals6 = vals[is.element(vals$time,re.for.corr[i,"date"]+5),]
      tmp.vals6 = tmp.vals6[which(tmp.vals6$bound=="mean"),]
      
      tmp.vals7 = vals[is.element(vals$time,re.for.corr[i,"date"]+6),]
      tmp.vals7 = tmp.vals7[which(tmp.vals7$bound=="mean"),]
      
      tmp.vals8 = vals[is.element(vals$time,re.for.corr[i,"date"]+7),]
      tmp.vals8 = tmp.vals8[which(tmp.vals8$bound=="mean"),]
      
      
      
      # take the average
      new.data = tmp.vals1
      new.re.for.corr.8days = re.for.corr[i,]
      for (j in seq(1,length(new.data)-2)){
        new.data[1,j] = (new.data[1,j] + tmp.vals2[1,j] + tmp.vals3[1,j] + tmp.vals4[1,j]  + tmp.vals5[1,j]  + tmp.vals6[1,j]   + tmp.vals7[1,j]   + tmp.vals8[1,j])/8
      }
      for (j in seq(1,length(re.for.corr)-1)){
        new.re.for.corr.8days[1,j] = (new.re.for.corr.8days[1,j] + re.for.corr[i+1,j] + re.for.corr[i+2,j] + re.for.corr[i+3,j])/4
      }
      if (i==1){
        weather.for.corr.8days = new.data
        re.for.corr.8days = new.re.for.corr.8days
      }else{
        weather.for.corr.8days = rbind(weather.for.corr.8days,new.data)
        re.for.corr.8days = rbind(re.for.corr.8days,new.re.for.corr.8days)
      }
    }
    
    
    # Average values over 10 days ---------------------------------
    
    
    remove(weather.for.corr.10days)
    remove(re.for.corr.10days)
    
    for (i in seq(1,length(re.for.corr$dat), 4) ){
      tmp.vals1 = vals[is.element(vals$time,re.for.corr[i,"date"]),]
      tmp.vals1 = tmp.vals1[which(tmp.vals1$bound=="mean"),]
      
      tmp.vals2 = vals[is.element(vals$time,re.for.corr[i,"date"]+1),]
      tmp.vals2 = tmp.vals2[which(tmp.vals2$bound=="mean"),]
      
      tmp.vals3 = vals[is.element(vals$time,re.for.corr[i,"date"]+2),]
      tmp.vals3 = tmp.vals3[which(tmp.vals3$bound=="mean"),]
      
      tmp.vals4 = vals[is.element(vals$time,re.for.corr[i,"date"]+3),]
      tmp.vals4 = tmp.vals4[which(tmp.vals4$bound=="mean"),]
      
      tmp.vals5 = vals[is.element(vals$time,re.for.corr[i,"date"]+4),]
      tmp.vals5 = tmp.vals5[which(tmp.vals5$bound=="mean"),]
      
      tmp.vals6 = vals[is.element(vals$time,re.for.corr[i,"date"]+5),]
      tmp.vals6 = tmp.vals6[which(tmp.vals6$bound=="mean"),]
      
      tmp.vals7 = vals[is.element(vals$time,re.for.corr[i,"date"]+6),]
      tmp.vals7 = tmp.vals7[which(tmp.vals7$bound=="mean"),]
      
      tmp.vals8 = vals[is.element(vals$time,re.for.corr[i,"date"]+7),]
      tmp.vals8 = tmp.vals8[which(tmp.vals8$bound=="mean"),]
      
      tmp.vals9 = vals[is.element(vals$time,re.for.corr[i,"date"]+8),]
      tmp.vals9 = tmp.vals9[which(tmp.vals9$bound=="mean"),]
      
      tmp.vals10 = vals[is.element(vals$time,re.for.corr[i,"date"]+9),]
      tmp.vals10 = tmp.vals10[which(tmp.vals10$bound=="mean"),]
      
      
      # take the average
      new.data = tmp.vals1
      new.re.for.corr.10days = re.for.corr[i,]
      for (j in seq(1,length(new.data)-2)){
        new.data[1,j] = (new.data[1,j] + tmp.vals2[1,j] + tmp.vals3[1,j] + tmp.vals4[1,j]  + tmp.vals5[1,j]  + tmp.vals6[1,j]   + tmp.vals7[1,j]   + tmp.vals8[1,j]    + tmp.vals9[1,j]    + tmp.vals10[1,j])/10
      }
      for (j in seq(1,length(re.for.corr)-1)){
        new.re.for.corr.8days[1,j] = (new.re.for.corr.8days[1,j] + re.for.corr[i+1,j] + re.for.corr[i+2,j] + re.for.corr[i+3,j]  + re.for.corr[i+4,j])/5
      }
      if (i==1){
        weather.for.corr.10days = new.data
        re.for.corr.10days = new.re.for.corr.10days
      }else{
        weather.for.corr.10days = rbind(weather.for.corr.10days,new.data)
        re.for.corr.10days = rbind(re.for.corr.10days,new.re.for.corr.10days)
      }
    }
    
    
    # Compute Correlations for temperature humidity and school days ---------------------------------
    
    
    
    weather.labels = labels(weather.for.corr.4days)
    plot_labels = c(6,7,24)
    
    library(Hmisc)
    
    for (pl in seq(1,length(plot_labels))){
      i = plot_labels[[pl]]
      weather.labels[[2]][i] 
      cor.vals = rcorr(re.for.corr$mean, weather.for.corr[,weather.labels[[2]][i]], type="pearson")
      new.corr.vals = data.frame(label=weather.labels[[2]][i], correlation=cor.vals$r["x","y"], p_value=cor.vals$P["x","y"], days=2, c=c, sub=sub)
      if (i==6 && c==1 && sub==1){
        corr.vals.days = new.corr.vals
      }else{
        corr.vals.days = rbind(corr.vals.days, new.corr.vals)
      }
    }
    
    for (pl in seq(1,length(plot_labels))){
      i = plot_labels[[pl]]
      weather.labels[[2]][i] 
      cor.vals = rcorr(re.for.corr.4days$mean, weather.for.corr.4days[,weather.labels[[2]][i]], type="pearson")
      new.corr.vals = data.frame(label=weather.labels[[2]][i], correlation=cor.vals$r["x","y"], p_value=cor.vals$P["x","y"], days=4, c=c, sub=sub)
      corr.vals.days = rbind(corr.vals.days, new.corr.vals)
    }
    
    for (pl in seq(1,length(plot_labels))){
      i = plot_labels[[pl]]
      weather.labels[[2]][i] 
      cor.vals = rcorr(re.for.corr.6days$mean, weather.for.corr.6days[,weather.labels[[2]][i]], type="pearson")
      new.corr.vals = data.frame(label=weather.labels[[2]][i], correlation=cor.vals$r["x","y"], p_value=cor.vals$P["x","y"], days=6, c=c, sub=sub)
      corr.vals.days = rbind(corr.vals.days, new.corr.vals)
    }
    
    for (pl in seq(1,length(plot_labels))){
      i = plot_labels[[pl]]
      weather.labels[[2]][i] 
      cor.vals = rcorr(re.for.corr.8days$mean, weather.for.corr.8days[,weather.labels[[2]][i]], type="pearson")
      new.corr.vals = data.frame(label=weather.labels[[2]][i], correlation=cor.vals$r["x","y"], p_value=cor.vals$P["x","y"], days=8, c=c, sub=sub)
      corr.vals.days = rbind(corr.vals.days, new.corr.vals)
    }
    
    
    for (pl in seq(1,length(plot_labels))){
      i = plot_labels[[pl]]
      weather.labels[[2]][i] 
      cor.vals = rcorr(re.for.corr.10days$mean, weather.for.corr.10days[,weather.labels[[2]][i]], type="pearson")
      new.corr.vals = data.frame(label=weather.labels[[2]][i], correlation=cor.vals$r["x","y"], p_value=cor.vals$P["x","y"], days=10, c=c, sub=sub)
      corr.vals.days = rbind(corr.vals.days, new.corr.vals)
    }
    
    
  }
  
  
  
  
  
  
  # plot the sampling proprotion
  p.samp[[sub]] = ggplot(data = sampling, aes(x=sampling.proportion)) +
    geom_histogram(aes(y=..density..), binwidth=0.001, colour="black", fill="white") + 
    xlab("proportion of sampled individuals")  +
    xlim(c(0.025,0.055)) +
    theme_minimal()+
    geom_vline(aes(xintercept=median(sampling.proportion, na.rm=T)),   # Ignore NA values for mean
               color="red", linetype="dashed", size=1)
  
}
require(grid)
require(gridExtra)
plot.sampling = do.call("grid.arrange",c(p.samp, ncol=5))
ggsave(plot=plot.sampling,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Supplement/sampling_proportion_all.pdf",width=15, height=7)
ggsave(plot=plot.sampling,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Supplement/sampling_proportion_all.png",width=15, height=7)

plot.weather = do.call("grid.arrange",c(p_wea, ncol=4))
ggsave(plot=plot.weather,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Supplement/bdsky_weather_all.pdf",width=20, height=12)
ggsave(plot=plot.weather,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Supplement/bdsky_weather_all.png",width=20, height=12)



p.pval <- ggplot() +
  geom_line(data=corr.vals.days, aes(x=days, y=p_value, color=as.character(c), group=interaction(sub,c))) +
  scale_y_log10()+
  theme_minimal()+
  scale_color_OkabeIto(  breaks = c("1", "2"),
                         labels = c("2016-11-01 until 2017-03-01", "2016-12-01 until 2017-02-01")
  ) +
  facet_wrap(label~.)+
  xlab("number of days data was averaged over") +
  ylab("p-value") +   labs(color = "time interval considered\nfor correlation")
plot(p.pval)

p.corrcoeff <- ggplot() +
  geom_line(data=corr.vals.days, aes(x=days, y=correlation, color=as.character(c), group=interaction(sub,c))) +
  theme_minimal()+
  scale_color_OkabeIto(  breaks = c("1", "2"),
                         labels = c("2016-11-01 until 2017-03-01", "2016-12-01 until 2017-02-01")
  ) +
  facet_wrap(label~.) +
  xlab("number of days data was averaged over") +
  ylab("correlation coefficient")  +   labs(color = "time interval considered\nfor correlation")
plot(p.corrcoeff)

ggsave(plot=p.pval,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Supplement/Correlation/p_value_subset.pdf", width=7, height=3)
ggsave(plot=p.pval,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Supplement/Correlation/p_value_subset.png", width=7, height=3)
ggsave(plot=p.corrcoeff,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Supplement/Correlation/corrcoeff_subset.png", width=7, height=3)
ggsave(plot=p.corrcoeff,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Supplement/Correlation/corrcoeff_subset.pdf", width=7, height=3)

