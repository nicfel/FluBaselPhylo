######################################################
######################################################
# combine the gene tree runs and run the mcc trees
######################################################
######################################################
library(ggplot2)
# needed to calculate ESS values
library(coda)
library(XML)
library("methods")


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
trees <- list.files(path="./bdskynucdiff/out/", pattern="rep0.*\\.trees$", full.names = TRUE)

system("rm -r bdskynucdiff/mcc")
system("mkdir bdskynucdiff/mcc")
system("rm -r bdskynucdiff/combined")
system("mkdir bdskynucdiff/combined")

# run log combiner
for (i in seq(1,length(trees))){
  in_command <- " -trees -resample 1000000 -burnin 1000000"
  for (j in seq(0,2)){
    in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), trees[i]), sep="")
  }

  out_command = gsub("rep0_", "", trees[i])
  out_command = gsub("out", "mcc", out_command)
  
  combined_out = gsub("mcc", "combined", out_command)
  
  
  system(paste("/Applications/beast1/bin/logcombiner ", in_command, combined_out, sep=" "))
  system(paste("/Applications/BEAST\\ 2.5.0/bin/treeannotator -burnin 0 -heights median ",combined_out, out_command, sep=" "))
}


logs <- list.files(path="./bdskynucdiff/out/", pattern="rep0.*\\.log$", full.names = TRUE)

# run log combiner
for (i in seq(1,length(logs))){
  system("rm tmp.trees")
  in_command <- "-resample 100000 -burnin 1000000"
  for (j in seq(0,2)){
    in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), logs[i]), sep="")
  }
  
  out_command = gsub("rep0_", "", logs[i])
  out_command = gsub("out", "combined", out_command)
  system(paste("/Applications/beast1/bin/logcombiner ", in_command, " ", out_command, sep=""))
}


