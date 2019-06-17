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
log <- list.files(path="./bdsky/out", pattern="*bdsky*.log", full.names = TRUE)

for (i in seq(1, length(log))){
  # Read in the MASCOT *.logs
  t <- read.table(log[i], header=TRUE, sep="\t")
  # take a 10 % burn in
  t <- t[-seq(1,ceiling(length(t$Sample)/10)), ]
  # get all labels
  all_labels = labels(t)
  new.rate <- data.frame(s=rep(i, length(t$Sample)))
  # calculate the R0's in real space at all time points
  for (j in seq(1, length(all_labels[[2]]))){
    if (startsWith(all_labels[[2]][[j]],"absolute")){
      absoluteR0 = t[,all_labels[[2]][[j]]]
    }
    if (startsWith(all_labels[[2]][[j]],"log")){
      logeR0 = t[,all_labels[[2]][[j]]]
      new.rate[, all_labels[[2]][[j]]] <- exp(logeR0)*absoluteR0
    }
  }
  
  if (i==1){
    rate=new.rate
  }else{
    rate = rbind(rate,new.rate)
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

# library(ggtree)
# require("OutbreakTools")
# tr <- read.annotated.nexus(file="all_trees.trees")
# 
# for (i in seq(1,length(tr$tip.label))){
#   # find correspondin edge
#   edge = which(tr$edge[,2]==i)
#   
#   new.labels <- data.frame(node=i, color=tr$annotations[[edge[1]]]$known)
#   if (i==1){
#     labels = new.labels
#   }else{
#     labels <- rbind(labels, new.labels)
#   }
# }
# 
# 
# for (i in seq(length(tr$tip.label)+1,length(tr$annotations))){
#   # find correspondin edge
#   edge = which(tr$edge[,1]==i)
#   
#   new.labels <- data.frame(node=i, color=tr$annotations[[edge[1]]]$known)
#   labels <- rbind(labels, new.labels)
# 
# }

# add age information
trait <- read.table("quest_labels.csv", sep=",", stringsAsFactor=F)
age <- read.csv("age_labels.csv", sep=",")

colnames(trait) <- sub("\\.$", "", colnames(trait))


color_vals =  colorRampPalette(c("red", "green"))(max(age$age)+1)

for (i in seq(1, length(age$taxa))){
  new.taxa <- data.frame(taxa=age[i,"taxa"], age=age[i,"age"], color = color_vals[age[i,"age"]+1])
  if (i==1){
    taxa = new.taxa
  }else{
    taxa = rbind(taxa, new.taxa)
  }
}


require(ggplot2)
require(ggtree)
x <- read.beast("all_trees.trees")
cols=scale_color(x, by="known", high="#000000", low = "#FFFFFF")
p_tree <- ggtree(x, right=TRUE, mrsd="2017-04-05",  color=cols, aes(size=known)) + theme_tree2() +
  geom_text(aes(x=max(x), label=NA), size=1, color=cols, hjust=-.3) +
  scale_x_continuous(breaks=c(2016.915301, 2017 , 2017.084932, 2017.161644), label = c("2016-12-01", "2017-01-01", "2017-02-01", "2017-03-01")) +
  # geom_segment(aes(xend=max(x)+.20, yend=y), linetype="dotted", size=.1, color=cols) +
  # theme(panel.grid.major   = element_line(color="black", size=.2),
  #       panel.grid.minor   = element_line(color="grey", size=.2),
  #       panel.grid.major.y = element_blank(),
  #       panel.grid.minor.y = element_blank()) + 
  scale_size(range=c(0,0.1) )

p_tree <- (p_tree) %<+% taxa + geom_tippoint(aes(color=age), shape=16, size=1.8) + scale_colour_gradientn(colours=c("red", "yellow","blue","green"))
gheatmap(p_tree, trait, offset=0, width=0.3, font.size=3, colnames_angle=90, hjust=1) 

ggsave(plot=gheatmap(p_tree, trait, offset=0, width=0.1, font.size=3, colnames_angle=90, hjust=1),"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/bdsky/tree.pdf",width=8, height=10)


# put the time axis in the right time


dat$t = dat$t-max(dat$t)
dat$t = dat$t*1.5/365
dat$t = dat$t + 2017.257534


p <- ggplot(data=dat) + geom_line(aes(x=t, y=mean)) +
  geom_line(aes(x=t, y=upper)) +
  geom_line(aes(x=t, y=lower)) +
  scale_x_continuous(breaks=c(2016.915301, 2017 , 2017.084932, 2017.161644), label = c("2016-12-01", "2017-01-01", "2017-02-01", "2017-03-01")) +
  # geom_segment(aes(xend=max(x)+.20, yend=y), linetype="dotted", size=.1, color=cols) +
  # theme(panel.grid.major   = element_line(color="black", size=.2),
  #       panel.grid.minor   = element_line(color="grey", size=.2),
  #       panel.grid.major.y = element_blank(),
  #       panel.grid.minor.y = element_blank())  
  theme_minimal() + xlab("time") + ylab("R0") + ylim(c(0,3))

plot(p)

ggsave(plot=p,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Global/bdsky.pdf",width=8, height=5)

