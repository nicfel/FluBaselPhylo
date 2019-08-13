# LocalClusters

For each of the initial clusters that were done based on nucleotide differences containing sequences from this study and from [gisaid.org](gisaid.org), a BEAST2 xml is made using a constant coalescent prior and the evolutionary rates estimated in the evolutionary rates folder.
For each logged iteration on the BEAST2 run, ancestral state reconstruction is then performed in order to assess which sequences are in the same local cluster, which is then used to dermine the number of introductions and the effective reproduction number over time.
