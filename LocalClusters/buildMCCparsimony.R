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



trees <- list.files(path="./typed/combined/", pattern="*\\.trees$", full.names = TRUE)

system("rm -r typed/mcc")
system("mkdir typed/mcc")

# run log combiner
for (i in seq(1,length(trees))){

  out_command = gsub("combined", "mcc", trees[i])
  
  combined_command = gsub(".trees",".combined.trees", out_command)
  combined_command = gsub("mcc//","combined/", combined_command)


  system(paste("/Applications/BEAST\\ 2.5.0/bin/treeannotator -burnin 0 -heights median", trees[i], out_command, sep=" "))
}

