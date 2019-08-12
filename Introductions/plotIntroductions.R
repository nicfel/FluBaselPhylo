library(ggplot2)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

t = read.table(file="./out/nr_introductions.csv", header=T, sep="\t")

p <- ggplot(t, aes(x=nr_introductions)) +
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") + xlab("number of sampled introductions")  +
  theme_minimal()+
  geom_vline(aes(xintercept=mean(nr_introductions, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)
plot(p)

ggsave(plot=p,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Clusters/introductions.pdf",width=5, height=5)

t = read.table(file="./out/nr_introductions_subset.csv", header=T, sep=",")

p <- ggplot(data=t) +
  geom_line(aes(x=subset,y=lower), binwidth=1, colour="grey", fill="white") + 
  geom_line(aes(x=subset,y=mean), binwidth=1, colour="black", fill="white") + 
  geom_line(aes(x=subset,y=upper), binwidth=1, colour="grey", fill="white") + 
  xlab("number of samples")  +
  ylab("number of sampled introductions")  +
  theme_minimal()
plot(p)

ggsave(plot=p,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Clusters/introductions_subset.pdf",width=3, height=3)

