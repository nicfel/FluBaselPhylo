clear 
% define the segments

% define the number of samples
use_nr_samples = 1;

segments = {'HA' 'M' 'NA' 'NP' 'NS' 'PA' 'PB1' 'PB2'};
cluster_names = cell(0,0);
% get the different clusters names
fasta = fastaread(['../Clusters/clusters/' segments{1} '.fasta']);
for j = 1 : length(fasta)
    tmp = strsplit(fasta(j).Header, '_');
    cluster_names{end+1} = tmp{1};
end
cluster_names = unique(cluster_names);
for j = 1 : length(segments)
    counter.(segments{j}) = 1;
end

for i = 1 : length(segments)
    fasta = fastaread(['../Clusters/clusters/' segments{i} '.fasta']);
    for j = 1 : length(fasta)
        tmp = strsplit(fasta(j).Header, '_');
        Data.(segments{i})(counter.(segments{i})) = fasta(j);
        counter.(segments{i}) = counter.(segments{i})+1;
    end
end

% check if everything's in order
name = cell(0,0);
for j = 1 : length(Data.(segments{1}))
    for k = 1 : length(segments)
        tmp = strsplit(Data.(segments{k})(j).Header,'|');
        name{end+1} = tmp{1};
    end
end    
if length(unique(name)) ~= length(Data.(segments{1}))
    error(unique(name))
end
disp('everything in order')

%

% check which sequences are not full genome
name = cell(0,0);
for j = 1 : length(segments)
    for i = 1 : length(Data.(segments{j}))    
        if isempty(strfind(Data.(segments{j})(i).Sequence, 'A'))
            tmp = strsplit(Data.(segments{j})(i).Header,'|');
            name{end+1} = tmp{1};
        end
    end
end
name = unique(name);

% delete all sequence with header in name
for j = 1 : length(segments)
    for i = length(Data.(segments{j})) : -1 : 1
        tmp = strsplit(Data.(segments{j})(i).Header,'|');
        if ~isempty(find(ismember(name,tmp{1})))
            Data.(segments{j})(i) = [];
        end
    end
end
%
% check if everything's in order
name = cell(0,0);
for j = 1 : length(Data.(segments{1}))
    for k = 1 : length(segments)
        tmp = strsplit(Data.(segments{k})(j).Header,'|');
        name{end+1} = tmp{1};
    end
end    
if length(unique(name)) ~= length(Data.(segments{1}))
    error(unique(name))
end
disp('everything in order')


%
all_headers = cell(0,0);
for j = 1 : length(segments)
    for i = 1:length(Data.(segments{j}))
        tmp = strsplit(Data.(segments{j})(i).Header,'|');
        all_headers{end+1,1} = tmp{1};
    end
end
uni_headers = unique(all_headers);
for i = 1 : length(uni_headers)
    if length(find(ismember(all_headers, uni_headers{i}))) ~= 8
        disp(length(find(ismember(all_headers, uni_headers{i}))))
        disp(uni_headers{i})
    end
end

%
use_samples = sort(randsample(length(uni_headers), use_nr_samples));
for i = length(use_samples) : -1 : 1
    uni_headers(use_samples(i)) = [];
end


for j = 1 : length(segments)
    for i = length(Data.(segments{j})) : -1 : 1
        tmp = strsplit(Data.(segments{j})(i).Header,'|');
        if ~isempty(find(ismember(uni_headers,tmp{1})))
            Data.(segments{j})(i) = [];
        end
    end
end


% check if everything's in order
name = cell(0,0);
for j = 1 : length(Data.(segments{1}))
    for k = 1 : length(segments)
        tmp = strsplit(Data.(segments{k})(j).Header,'|');
        name{end+1} = tmp{1};
    end
end    
if length(unique(name)) ~= length(Data.(segments{1}))
    error(unique(name))
end
disp('everything in order before xml')

system('rm -r nonaligned');
system('mkdir nonaligned');


fludbsegments = {'HA' 'M1' 'NA' 'NP' 'NS1' 'PA' 'PB1' 'PB2'};

% add new york full genomes second ref
newyork = fastaread('otherSequences/NewYork.fasta');
for j = 1 : length(newyork) 
    tmp = strsplit(newyork(j).Header, '|');
    tmp2 = strsplit(tmp{3}, '/'); 
    if length(tmp2) ==3
        seg_ind = find(ismember(fludbsegments, tmp{2}));
        Data.(segments{seg_ind})(length(Data.(segments{seg_ind}))+1) = newyork(j);
        Data.(segments{seg_ind})(length(Data.(segments{seg_ind}))).Header = ...
            sprintf('%s|flu|?|%s-%s-%s|northamerica|usa|usa|usacanada|?|?|?|?',...
            tmp{1}, tmp2{3}, tmp2{1}, tmp2{2});
    end
end

% add new york full genomes second ref
europe = fastaread('otherSequences/Europe.fasta');
for j = 1 : length(europe) 
    tmp = strsplit(europe(j).Header, '|');
    tmp2 = strsplit(tmp{3}, '/'); 
    if length(tmp2) ==3
        seg_ind = find(ismember(fludbsegments, tmp{2}));
        Data.(segments{seg_ind})(length(Data.(segments{seg_ind}))+1) = europe(j);
        Data.(segments{seg_ind})(length(Data.(segments{seg_ind}))).Header = ...
            sprintf('%s|flu|?|%s-%s-%s|northamerica|usa|usa|usacanada|?|?|?|?',...
            tmp{1}, tmp2{3}, tmp2{1}, tmp2{2});
    end
end



% add California full genomes second ref
california = fastaread('otherSequences/California.fasta');
for j = 1 : length(california) 
    tmp = strsplit(california(j).Header, '|');
    tmp2 = strsplit(tmp{3}, '/'); 
    if length(tmp2) ==3
        seg_ind = find(ismember(fludbsegments, tmp{2}));
        Data.(segments{seg_ind})(length(Data.(segments{seg_ind}))+1) = california(j);
        Data.(segments{seg_ind})(length(Data.(segments{seg_ind}))).Header = ...
            sprintf('%s|flu|?|%s-%s-%s|europe|europe|europe|europe|?|?|?|?',...
            tmp{1}, tmp2{3}, tmp2{1}, tmp2{2});
    end
end

for j = 1 : length(segments) 
    fastawrite(['nonaligned/' segments{j} '.fasta'], Data.(segments{j}))
end


disp('do the alignment by hand');

