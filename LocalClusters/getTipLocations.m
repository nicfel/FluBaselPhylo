% build tip to location table
clear
% read in cluster fasta file
fasta = fastaread('../Clusters/clusters/HA.fasta');

% read in the reagions mapping
f = fopen('../NonSequenceData/geo_regions.tsv');
c = 1;
while ~feof(f)
    tmp = strsplit(fgets(f));
    r1{c} = lower(tmp{1});
    r2{c} = lower(tmp{2});
    r3{c} = lower(tmp{3});c = c + 1;
end

%%

f = fopen('tipLocations.csv', 'w');
fprintf(f, 'name,region1,region2,region3,region4\n');
for i = 1:length(fasta)
    tmp = strsplit(fasta(i).Header, '|');
    name = strsplit(tmp{1},'_');
    if length(name) > 2
        newname = [strrep(tmp{1}, [name{1} '_'],'') '.' name{1}];
    else
        newname = [name{2} '.' name{1}];
    end
    
    clear region4
    if strcmp(tmp{7}, 'basel')
        region4 = 'basel';
    else
        region4 = tmp{5};
    end
    ind = find(ismember(r1,tmp{6}));
    
    fprintf(f, '%s,%s,%s,%s,%s,%s\n', newname, tmp{5}, tmp{6}, tmp{7}, region4, r3{ind});
end
fclose(f);