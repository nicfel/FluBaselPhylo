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


# Read in the percentiles
age <- read.table("./out/age_distribution.csv", header = T, sep=",",check.names = F)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the percentiles as heatmaps
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

breaks <- unique(age$age_group)
b=breaks[seq(1,length(breaks),1)]


tic.labels = c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80-89","90-99")

p_w <- ggplot(data=age, aes(x=age_group)) +
  geom_histogram(binwidth=1, colour="black", fill="white",stat="count") +
  scale_x_discrete(breaks=b, labels=tic.labels)+

  ylab("number of patients") +
  xlab("age group") +
  theme(
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1))+ theme(legend.position="none")
plot(p_w)


ggsave(plot=p_w,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Mixing/Age_Distribution.pdf",width=5, height=5)

