######################################################
######################################################
# combine the gene tree runs and run the mcc trees
######################################################
######################################################
library(ape)
# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# read in the tip locations
# locations <- read.csv(file="./tipLocations.csv", header=TRUE, sep=",")



trees <- list.files(path="./constcoalnucdiff/out/", pattern="rep0.*\\.trees$", full.names = TRUE)

system("rm -r constcoalnucdiff/combined")
system("mkdir constcoalnucdiff/combined")

system("rm -r constcoalnucdiff/mcc")
system("mkdir constcoalnucdiff/mcc")

# run log combiner
for (i in seq(1,length(trees))){
  in_command <- " -b 10 -log"
  for (j in seq(0,9)){
    in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), trees[i]), sep="")
  }
  
  out_command = gsub("rep0_", "", trees[i])
  out_command = gsub("out", "mcc", out_command)
  
  combined_command = gsub(".trees",".combined.trees", out_command)
  combined_command = gsub("mcc//","combined/", combined_command)
  combined_command = paste(" -o ", combined_command, sep="")
  
  # combine the trees
  system(paste("/Applications/BEAST\\ 2.5.0/bin/logcombiner", in_command, combined_command, "", sep=" "))
  
}


logs <- list.files(path="./constcoalnucdiff/out/", pattern="rep0.*\\.log$", full.names = TRUE)

# run log combiner
for (i in seq(1,length(logs))){
  in_command <- " -b 10 -log"
  for (j in seq(0,9)){
    in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), logs[i]), sep="")
  }
  
  out_command = gsub("rep0_", "", logs[i])
  out_command = gsub("out", "mcc", out_command)
  out_command = paste(" -o ", out_command, sep="")
  
  system(paste("/Applications/BEAST\\ 2.5.0/bin/logcombiner ", in_command, " ", out_command, sep=""))
}

# run log combiner
for (i in seq(1,length(trees))){
  out_command = gsub("rep0_", "", trees[i])
  out_command = gsub("out", "mcc", out_command)
  
  combined_command = gsub(".trees",".combined.trees", out_command)
  combined_command = gsub("mcc//","combined/", combined_command)
  
  
  
  system(paste("/Applications/BEAST\\ 2.5.0/bin/treeannotator -burnin 0 -heights median", combined_command, out_command, sep=" "))
}

