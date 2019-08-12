
clear

cutoff = 1.0;

% read the HA raxml tree file that is converted such that it does not have
% branch information
f = fopen('../TimeTree/time_tree_nolabels.tree');
while ~feof(f)
    line = fgets(f);
    if length(line)>1000
        tree = line;
    end
end

ptree = phytreeread(tree);
%%
% get connectivity
con = getmatrix(ptree);
nodenames = get(ptree, 'nodenames');
isBasel = false(length(nodenames),1);
dist = pdist(ptree, 'SquareForm', true);

% get which nodes were sampled in basel
for i = 1 : length(nodenames)
    if ~isempty(strfind(nodenames{i},'CH/Basel'))
        isBasel(i) = true;
    end
end

baselnames = nodenames(isBasel);
baselDist = dist(isBasel,isBasel);
%%
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




%%
% get non_basel sequences in that cluster
otherseqs_dist = dist(isBasel,:);
otherseqs_dist(otherseqs_dist>cutoff) = NaN;
otherseqs_dist(~isnan(otherseqs_dist)) = 1;

u_cluster = unique(cluster);
for i = 1 : length(u_cluster)
    % get the indices of basel sequences in the cluster
    is_in_cluster = find(cluster==u_cluster(i));
    indices = [];
    for j = 1 : length(is_in_cluster)
        indices = [indices find(otherseqs_dist(is_in_cluster(j),:)==1)];
    end
    Cluster_members{i} = nodenames(unique(indices));
end

%%
% build an xml with the full genomes of the influenza sequences that infers
% individual transmission trees for each initial transmission cluster
% read in the consensus sequences of each segment
fastas = dir('../Sequences/tmp/*.fasta');
for i = 1 : length(fastas)
    if length(fastas(i).name) < 12
        tmp = strsplit(fastas(i).name, '.');
        dat1 = fastaread(['../Sequences/tmp/' fastas(i).name]);
        name2 = strrep(fastas(i).name, 'M','MP');
        dat2 = fastaread(['../Sequences/gisaid/gisaid_' name2]);
        dat = [dat1;dat2];
        seg.(tmp{1}) = dat;
    end
end
%%

% specify sequences that are not to be used
dont_use.HA = {''};
dont_use.M = {''};
dont_use.NA = {''};
dont_use.NP = {''};
dont_use.NS = {''};
dont_use.PA = {''};
dont_use.PB1 = {''};
dont_use.PB2 = {''};

clear Data
segments = fieldnames(seg);
for i = 1 : length(Cluster_members)
    member_names = Cluster_members{i};
%     if length(member_names) > 1
        c = 0;
        for j = 1 : length(member_names)
            use = false;
            % check if date is valid
            for l = 1 : length(seg.HA)
                if ~isempty(strfind(seg.HA(l).Header,member_names{j}))
                    name = seg.HA(l).Header;
                    tmp = strsplit(seg.HA(l).Header, '|');
                    tmp2 = strsplit(tmp{4}, '-');
                    if length(tmp2)==3
                        use = true;
                        c = c + 1;
                    end
                    break;
                end
            end
            if ~use
                disp(['sequence ' member_names{j} ' doesn''t have sampling day']);
            else
                for k = 1 : length(segments)
                    sequenced = false;
                    for l = 1 : length(seg.(segments{k}))
                        if ~isempty(strfind(seg.(segments{k})(l).Header,member_names{j}))
                            sequenced = true;
                            if length(seg.(segments{k})(l).Sequence)>5000
                                sequenced = false;
                            end
                            if ~isempty(find(ismember(dont_use.(segments{k}), member_names{j})))
                                sequenced = false;
                            end
                            break;
                        end
                    end
                    if sequenced
                        Data.(['c' num2str(i)]).(segments{k})(c) = seg.(segments{k})(l);
                    else
                        Data.(['c' num2str(i)]).(segments{k})(c).Header = name;
                        Data.(['c' num2str(i)]).(segments{k})(c).Sequence = 'N';
                   end
                end
            end
        end
%     end
end
cluster_names = fieldnames(Data);

%% Check everything is correct
for i = 1 : length(cluster_names)
    for j = 1 : length(Data.(cluster_names{i}).(segments{1}))
        for k = 2 : length(segments)
            n1 = strsplit(Data.(cluster_names{i}).(segments{k-1})(j).Header,'|');
            n2 = strsplit(Data.(cluster_names{i}).(segments{k})(j).Header,'|');
            if ~strcmp(n1{1}, n2{1})
                error(Data.(cluster_names{i}).(segments{k-1})(j).Header)
            end
        end
    end    
end

%% align all sequences
warning('off','all')
system('rm -r clusters');
system('mkdir clusters');
for j = 1 : length(segments)
    delete('tmp.fasta');
    delete('tmpout.fasta');
    fprintf('align %s\n', segments{j});
    data_new = Data.(cluster_names{1}).(segments{j});
    for k = 1 : length(data_new)
        data_new(k).Header = [cluster_names{1} '_' data_new(k).Header];
    end
    fastawrite('tmp.fasta',data_new);
    for i = 2 : length(cluster_names)
        data_new = Data.(cluster_names{i}).(segments{j});
        for k = 1 : length(data_new)
            data_new(k).Header = [cluster_names{i} '_' data_new(k).Header];
        end
        fastawrite('tmp.fasta',data_new);
    end
    delete(sprintf('%s.fasta', segments{j}))
    [~,~] = system(sprintf('../Software/muscle3.8.31_i86darwin64 -diags1 -maxiters 1 -in tmp.fasta -out clusters/%s.fasta', segments{j}));
    new_dat = fastaread(sprintf('clusters/%s.fasta', segments{j}));
    % rebuild the clusters
    for i = 1 : length(cluster_names)
        c = 1;
        clear dat
        for k = 1 : length(new_dat)
            tmp = strsplit(new_dat(k).Header,'_');
            if strcmp(tmp{1},cluster_names{i})
                dat(c) = new_dat(k);
                dat(c).Header = strrep(dat(c).Header,[cluster_names{i} '_'],'');
                c = c + 1;
            end
        end
        Data.(cluster_names{i}).(segments{j}) = dat;
    end
    delete('tmp.fasta');
end

%%
for i = 1 : length(cluster_names)
    name = cell(0,0);
    for j = 1 : length(Data.(cluster_names{i}).(segments{1}))
        for k = 1 : length(segments)
            tmp = strsplit(Data.(cluster_names{i}).(segments{k})(j).Header,'|');
            name{end+1} = tmp{1};
        end
    end    
    if length(unique(name)) ~= length(Data.(cluster_names{i}).(segments{1}))
        error(unique(name))
    end
end
