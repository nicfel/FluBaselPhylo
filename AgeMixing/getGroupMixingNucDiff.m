function [] = getGroupMixingNucDiff(file_val, cut)
% calculates the number of patients from different or the same age group
% that are considered to be associated with one another

% this option is to run the whole thing on the cluster (called euler)
if exist('input_vals_group_mixing_nucdiff.mat')>0
    true_file_val = file_val;
    load('input_vals_group_mixing_nucdiff.mat')
    file_val = true_file_val;
else
    clear
    HA = fastaread('../Clusters/clusterNucleotidesDist/HA.fasta');
    basel_id = cell(0,0);
    for i = 1 : length(HA)
        if ~isempty(strfind(HA(i).Header, 'Basel'))
            tmp1 = strsplit(HA(i).Header, '|');
            tmp2 = strsplit(tmp1{1}, '/');
            basel_id{end+1,1} = tmp2{end};
        end
    end

    
    % read in the metadata
    f = fopen('../NonSequenceData/Master_table_processed.csv');
    c = 1; 
    % skip first line (headers)
    fgets(f);
    while ~feof(f)
        % replace empty cells with point
        line = strsplit(fgets(f),',');
        ind = find(ismember(basel_id, line{2}));
        if strcmp(strtrim(line{29}), '0') && ~isempty(ind)
            MetaData.id{c,1} = line{2}; 
            MetaData.age(c,1) = str2double(line{12});  
            if MetaData.age(c,1)<7
                MetaData.group(c,1) = 1;
            elseif MetaData.age(c,1)<18
                MetaData.group(c,1) = 2;
            elseif MetaData.age(c,1)>65
                MetaData.group(c,1) = 6;
            else
                MetaData.group(c,1) = 3;
            end               
                
            c = c+1;    
        end
    end
    fclose(f);
    
    % read in the household status
    f = fopen('../Questionnairs/HouselholdMember.csv');

    % skip first line (headers)
    fgets(f);
    while ~feof(f)
        % replace empty cells with point
        line = strsplit(fgets(f),',');

        age(1) = str2double(line{2});
        age(2) = str2double(line{3});

        age_ind = find(ismember(MetaData.id, line{1}));

        if length(age_ind)>0
            if str2double(line{6})==1
                if age(1)>0 || age(2)>0 && MetaData.group(age_ind)==3
                    MetaData.group(age_ind) = 4;
                elseif MetaData.group(age_ind)==3
                    MetaData.group(age_ind) = 5;
                end
            end
        end

    end
    fclose(f);



    % read in the pairwise distances between any two samples
    f = fopen('../LocalClusters/localclustersnucdiff/pairwise_distances.csv','r');
    c = 1;
    while ~feof(f)
        line = strsplit(fgets(f),',');
        distances = str2double(line(3:end));
        med = quantile(distances,0.5);
        % get the id of the individual
        id{1} = line{1};
        id{2} = line{2};

        ind1 = find(ismember(MetaData.id,id{1}));
        ind2 = find(ismember(MetaData.id,id{2}));
        
        % if one of the indices is empty, they are a sequences from the
        % sampe patient
        if ~isempty(ind1) && ~isempty(ind2)
        
            ind1 = ind1(1);
            ind2 = ind2(1);

            value_all(c,1) = MetaData.group(ind1);
            value_all(c,2) = MetaData.group(ind2);
            value_all(c,3) = quantile(distances,0.025);
            value_all(c,4) = med;
            value_all(c,5) = quantile(distances,0.975);
            value_all(c,6) = ind1;
            value_all(c,7) = ind2;

            dist(c) = med;
            c = c+1;
        end
    end
    fclose(f);

    inv_map = [2 1 3 4 5 7 6];
    value_dublicated = [value_all; value_all(:,inv_map)];

    % set upper and lower limits on whats considered to be associated
    c_off_upper = [0.05 0.1 0.15 0.2 0.3];
    c_off_lower = [0 0 0 0 0];

    % define the upper and lower bounds for age groups
    interval_sizes = [1];
    
    % define the number of permutations to use 
    nr_reps = 100000;

    save('input_vals_group_mixing_nucdiff')
    disp('only made input file')
    return
end
rng(file_val)

% override nr_reps
ints=1;
% set the cutoff levels    
upper = interval_sizes(ints)+1:1:7;
lower = upper - interval_sizes(ints);

% get all value_all entries that have distances between the lower and
% upper cutoff values
use_value_indices = find(value_dublicated(:,4)>=c_off_lower(cut) &...
                            value_dublicated(:,4)<c_off_upper(cut));

% only use the subset of pairs that have a distance withing the lower
% and upper bounds
value = value_dublicated(use_value_indices,:);

% make sanity tests
for a = 1 : length(upper)
    age_group_a = find(value(:,1)>=lower(a) & value(:,1) < upper(a));
    age_group_b = find(value(:,2)>=lower(a) & value(:,2) < upper(a));
    
    if length(age_group_a) ~= length(age_group_b)
        error('error in the number of pairs')
    end
    if length(unique(value(age_group_a,7))) ~= length(unique(value(age_group_b,6)))
        error('error in the number of pairs')
    end

    
end

clear nr_points_unique nr_points_pairs;
for a = 1 : length(upper)
    age_group_a = find(value(:,1)>=lower(a) & value(:,1) < upper(a));
    for b = 1 : length(upper)
        int = age_group_a(value(age_group_a,2)>=lower(b) & value(age_group_a,2) < upper(b));  
        % get all individuals involved
        vals = value(int,[6,7]);
        sorted_vals = sort(vals(:));
        sorted_vals(sorted_vals((1:end-1)) == sorted_vals((2:end)) ) = [];

        nr_points_unique(a,b) = length(sorted_vals);
        nr_points_pairs(a,b) = length(int);
    end
end



% sanity check if the matrices are symmetric
if ~issymmetric(nr_points_unique)
    error('unique not symmetric');
end
if ~issymmetric(nr_points_pairs)
    error('pairs not symmetric');
end

% get all ages as represented in the study
uni_age = MetaData.group;
[~,from,to] = unique(value(:,7));

% initialize the highs and lows 
high_vals = zeros(1,nr_reps);
low_vals = zeros(1,nr_reps);

vals_up = zeros(size(nr_points_unique));
vals_low = zeros(size(nr_points_unique));

vals_up_pairs = zeros(size(nr_points_unique));
vals_low_pairs = zeros(size(nr_points_unique));

% make a parallel loop over all permutations
for r = 1 : nr_reps 
    if mod(r,100)==0
        fprintf('%d %d %d\n', ints, cut, r)
    end

    % reshuffle age labels
    shuffel_age = uni_age(randsample(length(uni_age),length(uni_age)));
    
    age_group_a_prec = cell(length(upper),1);
    age_group_b_prec = cell(length(upper),1);
    P = cell(0,0);
    
    for a = 1 : length(upper)
        % get all individuals in the shuffeld age groups
        age_group_a = find(shuffel_age(value(:,6))>=lower(a) & shuffel_age(value(:,6)) < upper(a));
        age_partners = shuffel_age(value(age_group_a,7));
        for b = 1 : length(upper)
            % get the number of (shuffled) partners in the other age group
             int = age_group_a(age_partners>=lower(b) &....
                 age_partners < upper(b));              
             
%              % get the number of unique partners of all individuals
%              % involved
%              vals = value(int,[6,7]);
%              sorted_vals = sort(vals(:));
%              sorted_vals(sorted_vals((1:end-1)) == sorted_vals((2:end)) ) = [];
%              unis = length(sorted_vals);  
             
             pairs = length(int);

%              % check if the rand is higher or lower than the estimate
%              if unis < nr_points_unique(a,b)
%                  vals_up(a,b) = vals_up(a,b) + 1;
%              elseif unis > nr_points_unique(a,b)
%                  vals_low(a,b) = vals_low(a,b) + 1;
%              end        

              % check if the rand is higher or lower than the estimate
             if pairs < nr_points_pairs(a,b)
                 vals_up_pairs(a,b) = vals_up_pairs(a,b) + 1;
             elseif pairs > nr_points_pairs(a,b)
                 vals_low_pairs(a,b) = vals_low_pairs(a,b) + 1;
             end                     

        end
    end   
end


% % print to file
% f = fopen(sprintf('group_mixing_%d_%d_%d.csv',cut, interval_sizes(ints), file_val),'w');
% fprintf(f, '# lower = %f\n', c_off_lower(cut));
% fprintf(f, '# upper = %f\n', c_off_upper(cut));
% % make the header
% header = 'from,to,up,low';
% % print header
% fprintf(f,'%s\n',header);
% 
% % print percentiles
% for a = 1 : length(upper)
%     val_from = sprintf('%03d_%03d',lower(a),upper(a));
%     for b = a : length(upper)
%         val_to = sprintf('%03d_%03d',lower(b),upper(b));
%         fprintf(f,'%s,%s,%s,%s\n',val_from,val_to,num2str(vals_up(a,b)),num2str(vals_low(a,b)));
%     end
% end
% fclose(f);

% print the paris to file
f = fopen(sprintf('group_mixing_pairs_%d_%d_%d.csv',cut, interval_sizes(ints), file_val),'w');
fprintf(f, '# lower = %f\n', c_off_lower(cut));
fprintf(f, '# upper = %f\n', c_off_upper(cut));
% make the header
header = 'from,to,up,low';
% print header
fprintf(f,'%s\n',header);

% print percentiles
for a = 1 : length(upper)
    val_from = sprintf('%03d_%03d',lower(a),upper(a));
    for b = 1 : length(upper)
        val_to = sprintf('%03d_%03d',lower(b),upper(b));
        fprintf(f,'%s,%s,%s,%s\n',val_from,val_to,num2str(vals_up_pairs(a,b)),num2str(vals_low_pairs(a,b)));
    end
end
fclose(f);
        

end