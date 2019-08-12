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

# read in the cutoff values
filevals <- read.table("./out/file_values_group_mixing_pairs.csv", header = T, sep=",",check.names = F)

tic.labels = c("pre-school", "school","adults unknown","adults with kids","adults without kids","elderly")

# set and upper or lower limit for the association score (just for plotting reasons)
lim=2

# do the same plotting in a grid style tough
int_size <- unique(filevals$interval_size)
cut_offs <- unique(filevals$upper)

p <- list()

for (b in seq(1,length(cut_offs))){
  ind <- intersect(which(filevals$upper==cut_offs[b]),which(filevals$interval_size==int_size[1]))
  
  t = read.table(file=as.character(filevals[ind,"filename"]), header = T, sep=",",check.names = F)
  
  # t$percentile[which(t$percentile < 1.3010 & t$percentile > -0.5)] <- 0

  # set everything above the limiit to be the limit (due to color gradients)
  t$percentile[which(t$percentile< -lim)] <- -lim
  t$percentile[which(t$percentile> lim)] <- lim
  
  breaks = unique(t$from)
  ylabelstext=breaks[seq(1,length(breaks),1)]
  reverse_order = c("006_007", "005_006","004_005","003_004","002_003","001_002")
  
  # remove a part of the matrix
  for (i in seq(1,length(ylabelstext))){
    for (j in seq(1,length(ylabelstext))){
      if (j>i){
        t[which(t$from==ylabelstext[[i]] & t$to==ylabelstext[[j]]),]$percentile = NA
      }
    }
  }
  
  
  
  p[[b]] <- ggplot() +
    geom_tile(data=t, aes(ordered(to, levels=ylabelstext), ordered(from, levels=reverse_order), fill = percentile)) +
    scale_fill_gradientn(na.value="white", colors = c("#88419d","#b3cde3", rgb(1,1,1),  "#fdcc8a","#d7301f"),
                         limits = c(-lim,lim), 
                         breaks = c(-2, -1, 0, 1, 2), 
                         name="Association",
                         label=c(expression(paste("p" <= 0.01, sep="")),
                                 "p=0.1",
                                 "",
                                 "p=0.1",
                                 expression(paste("p" <= 0.01, sep=""))
                         )) +
    scale_x_discrete(breaks=ylabelstext, labels=tic.labels)+
    scale_y_discrete(breaks=ylabelstext, labels=tic.labels)+
    ylab("") +
    xlab("") +
    theme_minimal()+
    theme(
      # axis.ticks.x=element_blank(),
      # axis.text.x=element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title=element_text(size=14,face="bold"),
      legend.position = "none") +ggtitle(paste("cutoff value =", cut_offs[b], "years" ))
}

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend.plot <-p[[5]] + theme(legend.position="right"
                             )
legend <- g_legend(legend.plot) 
p[[6]] = legend
require(grid)
require(gridExtra)
plot <- do.call("grid.arrange",c(p, ncol=3))
ggsave(plot=plot,paste("/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Supplement/Group_Mixing_pairs_all.pdf", sep=""),width=10, height=8)
ggsave(plot=plot,paste("/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Supplement/Group_Mixing_pairs_all.png", sep=""),width=10, height=8)


# make the plot for interval size 10 years and cutoff 0.1 yeasr seperately
ind <- intersect(which(filevals$upper==0.1),which(filevals$interval_size==1))
t = read.table(file=as.character(filevals[ind,"filename"]), header = T, sep=",",check.names = F)
# t$percentile[which(t$percentile < 1.3010 & t$percentile > -1.3010)] <- 0
# set everything above the limiit to be the limit (due to color gradients)
t$percentile[which(t$percentile< -lim)] <- -lim
t$percentile[which(t$percentile> lim)] <- lim

breaks = unique(t$from)
ylabelstext=breaks[seq(1,length(breaks),1)]

# remove a part of the matrix
for (i in seq(1,length(ylabelstext))){
  for (j in seq(1,length(ylabelstext))){
    if (j>i){
      t[which(t$from==ylabelstext[[i]] & t$to==ylabelstext[[j]]),]$percentile = NA
    }
  }
}

p <- ggplot() +
  geom_tile(data=t, aes(ordered(to, levels=ylabelstext), ordered(from, levels=reverse_order), fill = percentile)) +
  scale_fill_gradientn(na.value="white", colors = c("#88419d","#b3cde3", rgb(1,1,1),  "#fdcc8a","#d7301f"),
                       limits = c(-lim,lim), 
                       breaks = c(-2, -1, 0, 1, 2), 
                       name="Association",
                       label=c(expression(paste("p" <= 0.01, sep="")),
                               "p=0.1",
                               "",
                               "p=0.1",
                               expression(paste("p" <= 0.01, sep=""))
                       )) +
  scale_x_discrete(breaks=ylabelstext, labels=tic.labels)+
  scale_y_discrete(breaks=ylabelstext, labels=tic.labels)+
  ylab("") +
  xlab("") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title=element_text(size=14,face="bold"))
  # + theme(legend.position="none")
plot(p)
ggsave(plot=p,paste("/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Age/Group_Mixing_pairs.pdf", sep=""),width=6, height=5)
ggsave(plot=p,paste("/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Age/Group_Mixing_pairs.png", sep=""),width=6, height=5)
