# Characterising the epidemic spread of Influenza A/H3N2 within a city through phylogenetics


Infecting large portions of the global population, seasonal influenza is a major burden on societies around the globe. 
In order to improve our understanding of how influenza is transmitted on a city scale, we collected an extremely densely sampled set of influenza sequences alongside patient metadata. 
We find that repeated introductions into the city drove the 2016/2017 influenza season and that downstream transmission dynamics correlated positively with temperature. 
Zooming into the local transmission outbreaks suggests that the elderly were to a large extent infected within their own transmission network, while school children drove the spread within the remaining transmission network. 
These patterns will be valuable to plan interventions combating the spread of influenza within cities given that similar patterns are observed for other influenza seasons and cities.

## AgeMixing

Contains the files to recreate figure 3B that computes the association between age groups based on random permutations of the group label to isolate. 

## Clusters

Contains the scripts for the initial clustering of sequences from Basel and around the globe into initial clusters that are later further characterized.

## Connections

Contains the scripts to recreate figures 2C and 3A.

## EvolutionaryRates

Contains the files to infer the evolutionary rates of influenza A/H3N2 over several season. 
These mean rate estimates are then later used for analyses where the evolutionary rates were fixed.

## Introductions

Computes the number of introductions based on the which sets of sequences from Basel are inferred to be in the same local clusters.
This is done for each iteration of the MCMC in the *LocalClusters* folder.
Based on that, a 95 % HPD for the number of introductions is calculated. 
This is done for several times with random subsets of the sequences from Basel in order to compute a the relation between introductions and number of samples as shown in figure 1C.

## LocalClusters

Contains the scripts that take the local clusters, based on the initial clustering using nucleotide differences, and creates the BEAST2 xml's to infer distributions of trees for each of these local clusters.
Additionally, this folder contains the scripts to calculate the internal locations of all nodes in these distributions of local trees using parsimony.
Also contains the scripts to compute the pairwise distances between any two sequences from Basel that are in the same initial clusters.

## Questionnaires

Contains the pdf version of the questionnaires (in german) distributed to patients, as well as

## R0

Contains the Reff through time analysis and the R scripts to recreate plot 1B and 1D-1G

## Sequences

Folder to get from rawreads to consensus sequences, also includes the gisaid sequences for which the acknowledgement table is provided

## Software

Contains the compiled bdsky package used for the Reff analysis

## TimeTree

Contains the pipeline to make the timetree shown in Figure 1. The pipeline uses raxml and timetree to build a timed phylogenetic tree and follows the (nextstrain.org)[nextstrain.org] pipeline.
