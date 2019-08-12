######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(coda)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


first = T

# get the names run log files
log <- list.files(path="./bdskysubset/out", pattern="rep0.8*bdsky.*log", full.names = TRUE)
for (i in seq(1, length(log))){
  print(i)
  # read in all three logs
  t.1 = read.table(log[i], header=TRUE, sep="\t")
  t.2 = read.table(gsub("rep0", "rep1", log[i]), header=TRUE, sep="\t")
  t.3 = read.table(gsub("rep0", "rep2", log[i]), header=TRUE, sep="\t")
  # take a 10 % burn in
  t.1 <- t.1[-seq(1,ceiling(length(t.1$Sample)/10)), ]
  t.2 <- t.2[-seq(1,ceiling(length(t.2$Sample)/10)), ]
  t.3 <- t.3[-seq(1,ceiling(length(t.3$Sample)/10)), ]
  # combne the log files
  t = rbind(t.1,t.2,t.3)
  # compute ess values
  ess = effectiveSize(t)
  if (ess[5]>50){
    # get teh hpd
    hpd = HPDinterval(as.mcmc(t$samplingProportion.1))
    # get the number of samples
    tmp = strsplit(strsplit(log[i], split="nr")[[1]][[2]], split="_")[[1]][[1]]
    new.samp = data.frame(nr=as.numeric(tmp), mean=mean(t$samplingProportion.1), lower = hpd[1,"lower"], upper = hpd[1,"upper"])
    if (first){
      samp = new.samp
      first = F
    }else{
      samp = rbind(samp, new.samp)
    }
  }else{
    print("low ESS")
  }
}

p = ggplot(samp) + 
  geom_point(aes(x=nr, y=mean))

plot(p)

