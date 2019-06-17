%% build input fasta for augur to get the initial HA tree
delete('ConsensusConverted/h3n2_ha.fasta');
system('cp ConsensusConverted/Ancestral_HA.fasta ConsensusConverted/h3n2_ha.fasta')
Data = fastaread('ConsensusConverted/backgroud_fullgenome_HA.fasta');

Data(1).Sequence = upper(Data(1).Sequence);
fastawrite('ConsensusConverted/h3n2_ha.fasta',Data);

Data = fastaread('ConsensusConverted/HA.fasta');
fastawrite('ConsensusConverted/h3n2_ha.fasta',Data);
