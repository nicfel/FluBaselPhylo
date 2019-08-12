library(graphics)
library(ape)
library("gridExtra")
library(ggplot2)
# library(treeio)
library(ggtree)
# library("OutbreakTools")

# clear workspace
rm(list = ls())



# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

tr <- read.tree("time_tree_ape.tree")
tr <- ladderize(tr, right = TRUE)

edges <- tr$edge
node_labs <- tr$node.label
tip_labs <- tr$tip.label


edge_labs <- matrix(, nrow = length(edges[,1]))
edge_labs[1:length(tip_labs)] <- tip_labs
edge_labs[(length(tip_labs)+1):(length(edges[,1])+1)] <- node_labs[1:(length(node_labs))]
edge_col = matrix(1, nrow = length(edge_labs)-1)


colors = c("#511EA8",
"#4138C3",
"#3E59CF",
"#4379CD",
"#4D92BE",
"#5AA5A8",
"#6BB18E",
"#7FB975",
"#96BD5F",
"#AFBD4F",
"#C5B945",
"#D8AE3E",
rgb(red=0.45, green=1,blue=1),
"#E39B39",
"#E67D33",
"#E2572B",
"#DC2F24")

use=FALSE
for (i in seq(1, length(tr$tip.label))){
  tmp = strsplit(tr$tip.label[i],split="_")
  new.labels <- data.frame(node=i, clade=tmp[[1]][2])
    
  new.taxlabels <- data.frame(taxa=tr$tip.label[i], basel="this study")

  if (as.numeric(tmp[[1]][3])==1){
    if (use){
      taxlabels = rbind(taxlabels, new.taxlabels)
    }else{
      taxlabels = new.taxlabels
      use=T
    } 
  }
  if (i==1){
    labels = new.labels
  }else{
    labels = rbind(labels, new.labels)
  }
}
for (i in seq(1, length(tr$node.label))){
  tmp = strsplit(tr$node.label[i],split="_")
  new.labels <- data.frame(node=i+length(tr$node.label), clade=tmp[[1]][2])
  labels = rbind(labels, new.labels)
}

# in here but not sure if needed
rownames(taxlabels) <- NULL
rownames(labels) <- NULL



colors = c("3c3" ="#511EA8",
           "3c3.A"="#4138C3",
           "3c2.A"="#3E59CF",
           "A3"="#4379CD",
           "A2"="#4D92BE",
           "A1"="#5AA5A8",
           "A1a"="#6BB18E",
           "A1b"="#7FB975",
           "A1b/135N"="#96BD5F",
           "A1b/135K"="#AFBD4F",
           "A4"="#C5B945")


# root height in true time, has to be done by hand
root_height = 2017.969403-6.0111615
  
  
p <- ggtree(tr) %<+% taxlabels + geom_tippoint(aes(size=basel), shape=16, color="red")
p <- p %<+% labels + aes(color=I(clade)) + scale_color_manual(values=colors) + guides(size=FALSE) + theme(legend.position="right")

plot(p)

p <- p + theme_minimal() +
  theme(panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=0.3),
        panel.grid.major.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(breaks = c(2012.0-root_height,2013-root_height,2014-root_height, 2015-root_height, 2016-root_height, 2017-root_height, 2018-root_height),
                     labels = seq(2012,2018,1)) #limits=c(2013.5-root_height, 2017.11-root_height))
  # geom_rect(xmin=2013.5-root_height, xmax=2014.5-root_height, ymin=0, ymax=Inf, fill="white", color="white")
  
p
ggsave(plot=p,"/Users/nicmuell/Documents/github/BaselFlu-Text/Figures/Global/Global_Tree.pdf",width=10, height=6)


