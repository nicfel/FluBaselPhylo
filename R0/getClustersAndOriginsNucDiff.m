% get which samples should be clustered together and which should be the
% prior on the origin
clear;fclose('all');

% read in the metadata (contains information about which individuals are
% dublicates
f = fopen('../NonSequenceData/Master_table_processed.csv');
c = 1; 
% skip first line (headers)
fgets(f);
while ~feof(f);
    % replace empty cells with point
    line = strsplit(fgets(f),',');
    if strcmp(strtrim(line{29}), '0')
        MetaData.id{c,1} = line{2}; 
        MetaData.samePatient(c,1) = str2double(line{11});
        MetaData.age(c,1) = str2double(line{12});  
        MetaData.samplingTime{c,1} = str2double(strsplit(line{6}, '.'));
        c = c+1;    
    end
end
fclose(f);





% set the random number seed
rng(1)

for nr_rep = 1 : 10

    % read in the lcoal clustering data
    f = fopen('../LocalClusters/localclustersnucdiff/parsimony_cluster_membership.csv');
    fgets(f);
    c = 1;
    cl_sets = cell(0,0);
    cl_sets_dim = cell(0,0);

    g = fopen(['localClustersAssignment/local_cluster_sets_nucdiff_nrrep' num2str(nr_rep) '.txt'], 'w');

    counter = 1;

    while ~feof(f)
        line = strtrim(fgets(f));
        % next clusters
        if strcmp(line, '-------')
            counter = counter+1;
            disp(counter);

            % get all the different sets
            uni_set = unique(cl_sets_dim);
            set_prob = zeros(1,length(uni_set));
            set_size = zeros(1,length(uni_set));


            % check the probability of the
            for i = 1 : length(cl_sets)
                for j = 1 : length(cl_sets{i})
                    ind = find(ismember(uni_set,cl_sets{i}{j}));
                    set_prob(ind) = set_prob(ind)+1;
                end
            end

            % get probs
            set_prob = set_prob/(c-1);


            % for each iterations, compute the average set prob
            avg_set_prob = zeros(length(cl_sets),1);
            set_size = zeros(length(cl_sets),1);
            ind_set = cell(length(cl_sets),1);
            for i = 1 : length(cl_sets)
                ind_set{i} = zeros(length(cl_sets{i}),1);
                for j = 1 : length(cl_sets{i})
                    ind = find(ismember(uni_set,cl_sets{i}{j}));
                    avg_set_prob(i) = avg_set_prob(i) + set_prob(ind);
                    ind_set{i}(j) = ind;
                end
                avg_set_prob(i) = avg_set_prob(i)/length(cl_sets{i});
                set_size(i) = length(cl_sets{i});
            end

            [max_val,~] = max(avg_set_prob);
            % get all the sets with the same prob
            ind_max = find(avg_set_prob==max_val);

            % if there is more than one iteration witht the highest set prob,
            % get the a random iteration
            if length(ind_max)>1
                ind_it = randsample(ind_max,1);
            else
                ind_it = ind_max;
            end

            ind_it = randsample(length(clusters),1);

            % keeps track of which sets to use
            use_set = false(1,length(uni_set));
            use_set(ind_set{ind_it}) = true;


            for i = 1 : length(uni_set)
                if use_set(i)
                    % go through every iterations and check for this set the
                    % maximal "origin" heights
                    height = zeros(0,0);
                    for j = 1 : length(cl_sets)
                        for k = 1 : length(cl_sets{j})
                            if strcmp(uni_set{i}, cl_sets{j}{k})
                                height(end+1) = clusters(j).heights(k);
                            end
                        end
                    end               

                    line = [uni_set{i} 'rep'];
                    fprintf(g, '%s\n', strrep(line, ',rep',...
                        ['|' num2str(set_prob(i))...
                         '|' num2str(quantile(height,0.025))...
                         '|' num2str(quantile(height,0.5))...
                         '|' num2str(quantile(height,0.975))]));
                end
            end

            % check that every sequence is assigned to a local cluster
            seqs_all = zeros(0,0);
            for i = 1 : length(uni_set)
                tmp = strsplit(uni_set{i}, ',');
                for j = 1 : length(tmp)-1
                    seqs_all(end+1) = str2double(tmp{j});
                end            
            end

            % make the same for the used subsets
            seqs_sub = zeros(0,0);
            for i = 1 : length(uni_set)
                if use_set(i)
                    tmp = strsplit(uni_set{i}, ',');
                    for j = 1 : length(tmp)-1
                        seqs_sub(end+1) = str2double(tmp{j});
                    end                  
                end
            end

            seqs_all = unique(seqs_all);
            seqs_sub = unique(seqs_sub);
            if length(seqs_all)~=length(seqs_sub)
                disp('some sequences were not assigned to a local cluster');
                disp(length(seqs_all));
                disp(length(seqs_sub));
    %         else
    %             valid = true;
            end





    %         % check if two sets are overlapping
    %         have_overlap = false(length(uni_set),length(uni_set));
    %         for a = 1 : length(uni_set)
    %             set1 = strsplit(uni_set{a},',');
    %             set_size(a) = length(set1)-1;
    %             for b = a+1 : length(uni_set)                
    %                 set2 = strsplit(uni_set{b},',');
    %                 interset = intersect(set1(1:end-1), set2(1:end-1));
    %                 if ~isempty(interset)
    %                     have_overlap(a,b) = true;
    %                     have_overlap(b,a) = true;
    %                 end
    %             end
    %         end
    %         
    %         % save the set prob
    %         set_prob = set_prob/(c-1);
    %         unchanged_set_prob = set_prob;        
    % 
    %         %%
    %         valid = false;
    %         while ~valid
    %             set_prob = unchanged_set_prob;
    %             % keeps track of which sets to use
    %             use_set = false(1,length(uni_set));
    %             % rank each set by it's probability, set the first one to be used,
    %             % then set the second one, if it's not overlapping to be used and
    %             % so on
    %             no_non_overlap = true;
    %             for i = 1 : length(uni_set)
    %                 [max_val,ind] = max(set_prob);
    %                 % get all the sets with the same prob
    %                 ind_max = find(set_prob==max_val);
    %                 
    %                 if max_val < 0.2
    %                     disp(ind_max)
    %                 end
    %                 
    %                 % if there is more than one with the same value, take the
    %                 % largest set
    %                 if length(ind_max)>1
    %                     % check if any of the sets overlap with previously
    %                     % added sets
    % %                     overlaps_withused = have_overlap(ind_max,use_set);
    % %                     overlaps_withused = sum(overlaps_withused,2)>0;
    % %                     
    % %                     ind_max(overlaps_withused)=[];
    % %                     ind_max
    %                     
    %                     ind = randsample(ind_max,1);                    
    %                 end
    % 
    %                 % check if ind overlaps with any previous sets
    %                 overlaps = have_overlap(ind,use_set);
    % 
    %                 if set_prob(ind) ==-1
    %                     break;
    % %                     error('set prob should not be here')
    %                 end
    %                 if sum(overlaps)==0
    %                     % the set has no overlap with any previous sets, therefore
    %                     % use it
    %                     disp(set_prob(ind))
    %                     use_set(ind) = true;
    %                     set_prob(ind) = -1;
    %                     % also set every set that overlaps with this one to -1
    %                     set_prob(have_overlap(ind,:)) = -1;
    %                 else
    %                     set_prob(ind) = -1;
    %                 end
    % 
    %             end
    % 
    %             for i = 1 : length(uni_set)
    %                 if use_set(i)
    %                     % go through every iterations and check for this set the
    %                     % maximal "origin" heights
    %                     height = zeros(0,0);
    %                     for j = 1 : length(cl_sets)
    %                         for k = 1 : length(cl_sets{j})
    %                             if strcmp(uni_set{i}, cl_sets{j}{k})
    %                                 height(end+1) = clusters(j).heights(k);
    %                             end
    %                         end
    %                     end               
    % 
    %                     line = [uni_set{i} 'rep'];
    %                     fprintf(g, '%s\n', strrep(line, ',rep',...
    %                         ['|' num2str(unchanged_set_prob(i))...
    %                          '|' num2str(quantile(height,0.025))...
    %                          '|' num2str(quantile(height,0.5))...
    %                          '|' num2str(quantile(height,0.975))]));
    %                 end
    %             end
    % 
    %             % check that every sequence is assigned to a local cluster
    %             seqs_all = zeros(0,0);
    %             for i = 1 : length(uni_set)
    %                 tmp = strsplit(uni_set{i}, ',');
    %                 for j = 1 : length(tmp)-1
    %                     seqs_all(end+1) = str2double(tmp{j});
    %                 end            
    %             end
    % 
    %             % make the same for the used subsets
    %             seqs_sub = zeros(0,0);
    %             for i = 1 : length(uni_set)
    %                 if use_set(i)
    %                     tmp = strsplit(uni_set{i}, ',');
    %                     for j = 1 : length(tmp)-1
    %                         seqs_sub(end+1) = str2double(tmp{j});
    %                     end                  
    %                 end
    %             end
    % 
    %             seqs_all = unique(seqs_all);
    %             seqs_sub = unique(seqs_sub);
    %             if length(seqs_all)~=length(seqs_sub)
    %                 disp('some sequences were not assigned to a local cluster');
    %                 disp(length(seqs_all));
    %                 disp(length(seqs_sub));
    %             else
    %                 valid = true;
    %             end
    %         end
            disp('found a valid combination of sets');

            %%

            % use the sets of sequences such that the probability of the
            % average set probability is maximal
            cl_sets = cell(0,0);
            cl_sets_dim = cell(0,0);
            c = 1;
        else
%             disp(line)
            cls = strsplit(line, ',');
            if isempty(line)
                break;
            end
            cl_sets{c} = cell(0,0);
            nonempty_cl_counter = 1;
            for i = 1 : length(cls)
                % split info and id's            
                tmp = strsplit(cls{i},':');
                height_tmp(i) = str2double(tmp{2});
                % get the id's of the members
                members_tmp{i} = strsplit(tmp{1},'|');
                % for each member, check if it's a dublicate
                for j = length(members_tmp{i}): -1 : 1
                    ind = find(ismember(MetaData.id, members_tmp{i}{j}));
                    % if is dublicate, prune the elements
                    if isempty(ind)
                        members_tmp{i}(j) = [];
                    end
                end
                clear nr
                if length(members_tmp{i})>0
                    for j = 1 : length(members_tmp{i})
                        tmp = strsplit(members_tmp{i}{j},'/');
                        tmp = strsplit(tmp{end},'.');
                        nr(j) = str2double(tmp{1});
                    end
                    members_nr{nonempty_cl_counter} = nr;
                    cl_sets{c}{nonempty_cl_counter} = sprintf('%d,', sort(nr));
                    cl_sets_dim{end+1} = sprintf('%d,', sort(nr));
                    nonempty_cl_counter = nonempty_cl_counter+1;
                end
            end
            clusters(c).heights = height_tmp;
            clusters(c).members = members_nr;
            c = c+1;
        end
    end
    fclose('all');
end