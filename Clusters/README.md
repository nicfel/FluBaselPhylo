# Clusters

In this Folder, a rough clustering of sequences is perfomed. The clustering is based on average nucleotide difference between sequences.
If any two sequences from this study have an average nucleotide difference per site of less than 0.0025, then they are assigned to the same cluster.
Clustering is therefore performed such that the maximal "link" distance between any two sequences from this study is less than 0.0025 nucleotide differences per site.
Then, all other sequences from [gisaid.org](gisaid.org) are added if their distance to any of the sequences from this study is less than  0.0025 nucleotide differences per site.
The factor of 2 is simply added in order to get cluster with number of sequences that are manageable for the downstream analyses
