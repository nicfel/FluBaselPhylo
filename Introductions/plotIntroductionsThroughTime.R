library(ggplot2)
library(lubridate)
library(coda)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# read in the file that contains all the sampling time information
sampling.info.file = read.table("../NonSequenceData/Master_table_processed.csv",
                                     header=T, sep=",")

# get the names run log files
log <- list.files(path="../R0/bdskynucdiff/combined", pattern="*heights.*.log", full.names = TRUE)
first = T
# for (i in seq(1, length(log))){
for (i in seq(1, 1)){
    
  tmp = strsplit(log[[i]], split="subs")
  tmp = gsub(".log","", tmp[[1]][[2]])
  
  # get the members of each cluster
  cluster.assignment.file = read.table(paste("../R0/localClustersAssignment/local_cluster_sets_nucdiff_nrrep", tmp, ".txt" , sep=""),
                                       header=F, sep="|")
  
  # Read in the tree heights
  t <- read.table(log[i], header=TRUE, sep="\t")
  # get all labels
  all_labels = labels(t)
  # loop over all clusters
  for (j in seq(1,length(cluster.assignment.file$V1))){
    # get all the samples in the cluster
    tmp = strsplit(as.character(cluster.assignment.file[j,"V1"]), split=",")
    
    mrsi = decimal_date(as.Date("2000-01-01"))
    
    for (k in seq(1,length(tmp[[1]]))){
      sampling.time = strsplit(as.character(sampling.info.file[which(sampling.info.file$Labor.Nummer==tmp[[1]][[k]]), "Probeneingang"]), split="\\.")
      new.date = decimal_date(as.Date(paste(paste("20",sampling.time[[1]][[3]],sep=""),sampling.time[[1]][[2]],  sampling.time[[1]][[1]], sep="-")))
      
      if (new.date>mrsi){
        mrsi=new.date
      }
      
    }
    
    # if the number of members in the cluster is larger than 1, get the tree height
    if (length(tmp[[1]])>1){
      name = paste("TreeHeight", j, "HA", sep=".")
      height.distribution = mrsi-t[!is.na(t[,name]),name]
      hpd = HPDinterval(as.mcmc(height.distribution))
      new.intro = data.frame(mean=mean(height.distribution), upper=hpd[1,"upper"], lower= hpd[1,"lower"], runnr=i)
    }else{
      new.intro = data.frame(mean=mrsi, upper=mrsi, lower=mrsi, runnr = i)
    }
    
    if (first){
      intro = new.intro
      first =F
    }else{
      intro = rbind(intro, new.intro)
    }
  }
}

# sort the entries
indices = sort(intro$mean, index.return=TRUE)$ix

intro = intro[indices,]
intro$y = seq(1,length(intro$mean))

p <- ggplot() +
  geom_errorbarh(data=intro, aes(y=y, xmin=lower, xmax=upper))
plot(p)


p <- ggplot() +
  geom_histogram(data=intro, aes(mean), binwidth=3.5/365)
plot(p)