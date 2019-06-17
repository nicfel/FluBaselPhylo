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
log <- list.files(path="./bdskynucdiff/out", pattern="*heights.*subs1.log", full.names = TRUE)

for (i in seq(1, length(log))){
  # Read in the MASCOT *.logs
  t <- read.table(log[i], header=TRUE, sep="\t")
  # take a 10 % burn in
  t <- t[-seq(1,ceiling(length(t$Sample)/10)), ]
  # get all labels
  all_labels = labels(t)
  new.height <- data.frame(s=rep(i, length(t$Sample)))
  # calculate the R0's in real space at all time points
  for (j in seq(1, length(all_labels[[2]]))){
    if (startsWith(all_labels[[2]][[j]],"TreeHeight")){
      new.height[, all_labels[[2]][[j]]] <- t[,all_labels[[2]][[j]]]
    }
  }
  
  if (i==1){
    height=new.height
  }else{
    height = rbind(height,new.height)
  }
  
}

# Put Everything into 1 row
for (i in seq(2, length(height))){
  new.combined <- data.frame(height=mean(height[,i]*365))
  if (i==2){
    combined=new.combined
  }else{
    combined = rbind(combined, new.combined)
  }
}

p <- ggplot(combined, aes(x=height)) +
  geom_histogram(aes(y=..density..), binwidth=5, colour="black", fill="white") + xlab("tree heights of local clusters in days")  +
  theme_minimal()+
  geom_vline(aes(xintercept=median(height, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)
plot(p)

ggsave(plot=p,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Clusters/heights.pdf",width=5, height=5)


# read in the mrca and persistance estimates from simulations
persist <- read.table("../Persistence/persistence.csv", header=TRUE, sep=",")

  
p_persist <- ggplot() + 
  geom_density(data=persist, aes(x=mrca*365, y=..density..,colour="simulated MRCA", fill="simulated MRCA") , alpha=0.25) 

p_persist <- p_persist +
  geom_density(data=combined, aes(x=height, y=..density.., colour="inferred MRCA", fill="inferred MRCA"), alpha=0.25)  
p_persist <- p_persist +
  geom_density(data=persist, aes(x=persist*365, y=..density.., colour="simulated persistence", fill="simulated persistence"), alpha=0.5)+
  xlab("days") +
  theme_minimal()+
  geom_vline(data=persist, aes(xintercept=median(mrca*365, na.rm=T)),   # Ignore NA values for mean
             color=col4, linetype="dashed", size=2)+
  geom_vline(data=combined, aes(xintercept=median(height, na.rm=T)),   # Ignore NA values for mean
             color=col0, linetype="dashed", size=1)+
  geom_vline(data=persist, aes(xintercept=median(persist*365, na.rm=T)),   # Ignore NA values for mean
             color=col2, linetype="dashed", size=1) + 
  scale_color_manual(values=c("simulated MRCA"=col4,"inferred MRCA"=col0,"simulated persistence"=col2))+
  scale_fill_manual(values=c("simulated MRCA"=col4,"inferred MRCA"=col0,"simulated persistence"=col2))

plot(p_persist)

ggsave(plot=p_persist,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Clusters/persistance.pdf",width=7, height=5)
