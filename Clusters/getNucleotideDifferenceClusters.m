%% load the input data
clear
load('input_distance.mat')
% combined the csv file into one matrix
for i = 1 : length(basel_ids)
    f = fopen(['diffout/dist_' num2str(i) '.csv']);
    t = textscan(f, '%f');
    fclose(f);
    baselDist_tmp(i,:) = t{1}';
end

cutoff = 0.0025;

% get the indices of the basel sequences
basel_inds = zeros(0,0);
for i = 1 : length(basel_ids)
    basel_inds(i) = find(ismember(uni_id, basel_ids{i}));
end
    

baselDist = baselDist_tmp(:,basel_inds);
% make a cutoff
baselDist(baselDist>cutoff) = NaN;
baselDist(~isnan(baselDist)) = 1;

% get the sequences that are in the same initial clusters
cluster = zeros(size(baselDist,2),1);
clusternr = 1;
clustered = zeros(0,0);
for i = 1 : size(baselDist,1)
    ind = find(clustered==i);
    % sample is not found in any cluster thus far
    if isempty(ind)
        % get all samples connected to the current cluster
        new_clusters = [i, find(baselDist(i,:)==1)];                        
        has_new_individuals = true;
        
        while has_new_individuals
            has_new_individuals = false;
            for j = 1 : length(new_clusters)
                new_neighbours = find(baselDist(new_clusters(j),:)==1);
                new_new_cluster = unique([new_clusters new_neighbours]);
                if length(new_new_cluster) ~= length(new_clusters)
                    has_new_individuals = true;
                    new_clusters = new_new_cluster;
                    break;
                end
            end
        end
        
        clustered = [clustered, new_clusters];        
        cluster(new_clusters)=clusternr;
        clusternr = clusternr + 1;           
    end
end

%% add non_basel sequences in that cluster
otherseqs_dist = baselDist_tmp;
otherseqs_dist(otherseqs_dist>cutoff/2) = NaN;
otherseqs_dist(~isnan(otherseqs_dist)) = 1;

clear Cluster_members
u_cluster = unique(cluster);
for i = 1 : length(u_cluster)
    % get the indices of basel sequences in the cluster
    is_in_cluster = find(cluster==u_cluster(i));
    indices = [];
    
    for j = 1 : length(is_in_cluster)
        indices = [indices basel_inds(is_in_cluster(j)) find(otherseqs_dist(is_in_cluster(j),:)==1)];
    end
    ids_to_add = uni_id(unique(indices));
    % remove ids with insufficient time label
    for j = length(ids_to_add):-1:1
        tmp = strsplit(ids_to_add{j}, '|');
        tmp2 = strsplit(tmp{4}, '-');
        if length(tmp2)~=3
            disp('rem')
            ids_to_add(j) = [];
        elseif ~isempty(strfind(tmp{4}, 'unknown'))
            ids_to_add(j) = [];
        end
    end
    Cluster_members{i} = ids_to_add;
end

%%

% print clusters based on nucleotide difference to files
for i = 1 : length(Data)
    count = 1;
    clear FastaData
    for c = 1 : length(u_cluster)
        for k = 1 : length(Cluster_members{c})
            ind = find(ismember(id{i}, Cluster_members{c}{k}));
            if isempty(ind)
                addData = Data{i}(1);
                addData.Header = Cluster_members{c}{k};
                addData.Sequence = strrep(addData.Sequence,'A','N');
                addData.Sequence = strrep(addData.Sequence,'T','N');
                addData.Sequence = strrep(addData.Sequence,'G','N');
                addData.Sequence = strrep(addData.Sequence,'C','N');
            else
                addData = Data{i}(ind);
            end
            FastaData(count) = addData;
            FastaData(count).Header = ['c' num2str(u_cluster(c)) '_' FastaData(count).Header];
            count = count+1;
        end        
    end
    fastawrite(['clusterNucleotidesDist/' strrep(fastas(i).name, 'afasta','fasta')], FastaData)
end
