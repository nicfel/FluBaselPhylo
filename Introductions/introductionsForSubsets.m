% analyze the source locations of introductions into basel
clear

% read in the master table to avoid using more than one sample per patient
f = fopen('../NonSequenceData/Master_table_processed.csv');
c = 1; 
% skip first line (headers)
fgets(f);
while ~feof(f)
    % replace empty cells with point
    line = strsplit(fgets(f),',');
    if strcmp(strtrim(line{29}), '0')
        MetaData.id{c,1} = line{2}; 
        c = c+1;    
    end
end
fclose(f);


%% read in the analysis
f = fopen('../LocalClusters/localclustersnucdiff/parsimony_cluster_membership.csv');

% get the locations;
locations = strsplit(fgets(f),',');

max_size = 500;

% keeps track of the id's that are dublicates
remove = cell(0,0);

cluster_cluster_size = sparse(0,max_size);
cluster_nr_introductions = zeros(1,max_size);

cluster_sizes = cell(0,0);

c_glob = 1;
c = 1;
while ~feof(f);
    line = strtrim(fgets(f));
    % check if its the end of the file
    if strcmp(line, '-------')
        c = 1;
        cluster_size(c_glob,:) = full(mean(cluster_cluster_size));
        
        cluster_sizes{c_glob} = cell(0,0);
        for i = 1 : size(cluster_cluster_size,1)
            non_zero = find(cluster_cluster_size(i,:));
            
            sizes = zeros(0,0);
            for j = 1:length(non_zero)
                sizes = [sizes repmat(non_zero(j),1,cluster_cluster_size(i,non_zero(j)))];
            end            
            cluster_sizes{c_glob}{i} = sizes;
        end        
       
        cluster_cluster_size = sparse(0,max_size);
        
        nr_introductions(c_glob,:) = mean(cluster_nr_introductions);
        cluster_nr_introductions = zeros(1,max_size);
        
        
        
        c_glob = c_glob + 1;
    else
        cluster = strsplit(line, ',');
        cluster_cluster_size(c,:) = sparse(1,max_size);
        nr_intros = 0;
        for i = 1 : length(cluster)
            % split info and id's
            tmp1 = strsplit(cluster{i}, ':');
            tmp2 = strsplit(tmp1{1}, '|');
            % check if some are dublicates
            % form the same patient
            cl_size = 0;
            for j = 1 : length(tmp2)
                ind = find(ismember(MetaData.id, tmp2{j}));
                if ~isempty(ind)                    
                    cl_size = cl_size+1;
                end
            end
            if cl_size>0
                nr_intros = nr_intros+1;
                cluster_cluster_size(c,cl_size) = cluster_cluster_size(c,cl_size)+1;
            end
        end
        % Can happen in case of multiple infections
        if nr_intros>0
            cluster_nr_introductions(c,nr_intros) = 1;
        else
            cluster_nr_introductions(c,:) = zeros(1,max_size);
        end
        c = c+1;
    end    
end

fclose(f);



%% calculate the sum of introductions
rng(1234569);

probs = nr_introductions(1,:);
nr_reps = 10000;

for r = 1 : nr_reps
    clear nr
    for i = 1 : size(nr_introductions,1)
        if sum(nr_introductions(i,:))>0
            nr(i) = randsample(length(nr_introductions(i,:)),1,true,nr_introductions(i,:));
        else
            nr(i) = 0;
        end
    end
    intros(r) = sum(nr);
end

f = fopen('out/nr_introductions.csv', 'w');
fprintf(f, 'nr_introductions\n');
for i = 1 : length(intros)
    fprintf(f, '%d\n', intros(i));
end
fclose(f)


%%
% combine the iterations
all_cl_size = cell(0,0);
for j = 1 : length(cluster_sizes{1})
    tmp = zeros(0,0);
    for i = 1 : length(cluster_sizes)
        tmp = [tmp cluster_sizes{i}{j}];
    end   
    all_cl_size{j,1} = tmp;
end

% take subsets 
subset = 1:1:663;

nr_subset_intros = zeros(length(subset), length(all_cl_size));
for i = 1 : length(subset)
    disp(i)
    for j = 1 : length(all_cl_size)
        tmp = all_cl_size{j};
        for k = 1 : sum(tmp)-subset(i)
            ind = randsample(length(tmp),1,true,tmp);
            tmp(ind) = tmp(ind)-1;
        end     
        nr_subset_intros(i,j) = length(tmp(tmp>0));
    end
end
    
f = fopen('out/nr_introductions_subset.csv', 'w');
fprintf(f, 'subset,lower,mean,median,upper\n');
for i = 1 : length(subset)
    fprintf(f, '%d,%f,%f,%f,%f\n', subset(i),quantile(nr_subset_intros(i,:),0.025),mean(nr_subset_intros(i,:)),quantile(nr_subset_intros(i,:),0.5),quantile(nr_subset_intros(i,:),0.975));
end
fclose(f)

% calculate the probability of adding a new introduction

f = fopen('out/prob_new_introduction.csv', 'w');
fprintf(f, 'subset,mean\n');
for i = 2 : length(subset)
    fprintf(f, '%d,%f\n', subset(i),mean(nr_subset_intros(i,:))-mean(nr_subset_intros(i-1,:)));
end
fclose(f)


