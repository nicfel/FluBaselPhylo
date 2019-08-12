% converts the resolved time tree to something readible by ape
clear
% read in the time resolved tree
f = fopen('time_tree_resolved.tree'); clear name loc
trees = cell(0,0);
tree_name = cell(0,0);
phytrees = cell(0,0);
print_line = true;
while ~feof(f)
    full_line = fgets(f);
    line = strsplit(strtrim(full_line));
    if length(line)==2
        val = str2double(line{1});
        if ~isnan(val)
            tmp = strrep(line{2},'''','');
            name{val,1} = strrep(tmp,',','');
        end
    elseif length(line)>4
        tmp = strrep(full_line, ' ','');
        tmp = strsplit(tmp, '[&R]');        
        
        nr_nodes = strfind(tmp{end}, ')');
        trees = sprintf(strrep(tmp{end},')',')n%d'),1:length(nr_nodes));      
        trees = strrep(trees, ');',sprintf(')n%d;',length(nr_nodes)+1));
        
        leafvals1 = regexp(trees, '\,(\d*)\[\&(.*?)\]:', 'match');
        leafvals2 = regexp(trees, '\((\d*)\[\&(.*?)\]:', 'match');
        nodevals = regexp(trees, 'n(\d*)\[\&(.*?)\]:', 'match');
        sampling_times_tmp1 = regexp(trees, '\,(\d*)\[\&num_date="(\d*).(\d*)"', 'match');
        sampling_times_tmp2 = regexp(trees, '\((\d*)\[\&num_date="(\d*).(\d*)"', 'match');
        trees = regexprep(trees, '\[\&(.*?)\]:', ':');

    end  
end
fclose(f);

% read in the tree
ptree = phytreeread(trees);
% get the values of the samplign times
sampling_times_tmp = [sampling_times_tmp1 sampling_times_tmp2];
sampling_times = zeros(length(sampling_times_tmp),1);
for i = 1 : length(sampling_times_tmp)
    tmp1 = strsplit(sampling_times_tmp{i}, '[');
    tmp1 = strrep(tmp1, ',', '');
    tmp1 = strrep(tmp1, '(', '');
    tmp2 = strsplit(sampling_times_tmp{i}, '"');
    sampling_times(str2double(tmp1{1})) = str2double(tmp2{2});
    
end


%% prune the outgroups
leafnames = get(ptree, 'leafnames');

prune_leafes = false(length(leafnames),1);
basel_leafes = false(length(leafnames),1);

sampling_times(sampling_times==0) = [];
count = 0;
for i = 1 : length(name)
    tmp = strsplit(name{i}, '/');
    
    if ~isempty(strfind(name{i}, 'A/CH/Basel'))
        ind = find(ismember(leafnames, num2str(i)));
        basel_leafes(ind) = true;
    else
        ind = find(ismember(leafnames, num2str(i)));
        
        st = sampling_times(str2double(leafnames{ind}));
        if st<2015.5 || st>2018.5
            prune_leafes(ind) = true;
        end
%         if binornd(1,0.8)
%             prune_leafes(ind) = true;
% %             sampling_times(str2double(leafnames{ind})) = 0;
%         end
    end
    
end

sum(prune_leafes==0)
%%
% also prunes two Ohio leafes that have ancestry way back
prune_leafes(1) = true;
prune_leafes(20) = true;

pruned_tree = prune(ptree, find(prune_leafes));

% get the newick string back
newick = getnewickstr(pruned_tree, 'BranchNames', true);

%% add the leave and node information
clear replace_nodes1 replace_nodes2
for i = 1 : length(nodevals)
    tmp = strsplit(nodevals{i}, '[');
    tmp = strrep(tmp{1}, ',', '');
    
    splits = strsplit(nodevals{i}, '"');    
    clade = splits{end-1};
    replace_nodes1{i} = tmp;
    replace_nodes2{i} = clade;
end
replace_leaves = cell(0,0);
for i = 1 : length(leafvals1)
    tmp = strsplit(leafvals1{i}, '[');
    tmp = strrep(tmp{1}, ',', '');    
    splits = strsplit(leafvals1{i}, '"');    
    clade = splits{end-1};
    replace_leaves{end+1,1} = tmp;
    replace_leaves{end,2} = clade;
end
for i = 1 : length(leafvals2)
    tmp = strsplit(leafvals2{i}, '[');
    tmp = strrep(tmp{1}, '(', '');    
    splits = strsplit(leafvals2{i}, '"');    
    clade = splits{end-1};
    replace_leaves{end+1,1} = tmp;
    replace_leaves{end,2} = clade;
end

% get the unique clades
unique_clades = unique([unique(replace_leaves(:,2)); unique(replace_nodes2(:))]);

%% add the labels 
print_tree = newick;
print_tree_sanity = newick;

for i = 1:length(replace_nodes1)
    if mod(i,1000)==1
        disp(i)
    end
    clade_nr = find(ismember(unique_clades, replace_nodes2{i}));
    print_tree = strrep(print_tree, [')' replace_nodes1{i} ':'],...
        [')' replace_nodes1{i} '_' replace_nodes2{i} '_0:']);
    print_tree_sanity = strrep(print_tree_sanity, [')' replace_nodes1{i} ':'],...
        [')[&clade=' replace_nodes2{i} ']:']);
end

use_leaf_names = cell(0,0);

for i = 1:length(replace_leaves)
    if mod(i,1000)==1
        disp(i)
    end
    val = str2double(replace_leaves{i,1});
    clade_nr = find(ismember(unique_clades, replace_leaves{i,2}));

    
    print_tree = strrep(print_tree, [',' replace_leaves{i,1} ':'],...
        [',' replace_leaves{i,1} '_' replace_leaves{i,2} '_' num2str(basel_leafes(val)) ':']);
    print_tree = strrep(print_tree, ['(' replace_leaves{i,1} ':'],...
        ['(' replace_leaves{i,1} '_' replace_leaves{i,2} '_' num2str(basel_leafes(val)) ':']);
    print_tree_sanity = strrep(print_tree_sanity, [',' replace_leaves{i,1} ':'],...
        [',' replace_leaves{i,1} '[&clade=' replace_leaves{i,2} ',basel=' num2str(basel_leafes(val)) ']:']);
    print_tree_sanity = strrep(print_tree_sanity, ['(' replace_leaves{i,1} ':'],...
        ['(' replace_leaves{i,1} '[&clade=' replace_leaves{i,2} ',basel=' num2str(basel_leafes(val)) ']:']);
end


%%


f = fopen('time_tree_ape.tree', 'w');
fprintf(f, 'tree TREE1 = [&R] %s\n', regexprep(print_tree, ')n(\d*);',');'));
fclose(f);

%%

leafnames = get(pruned_tree, 'leafnames');
leafnames_nr = str2double(leafnames);

f = fopen('time_tree_sanity.tree', 'w');
fprintf(f, '#NEXUS\n');
fprintf(f, 'Begin taxa;\n');
fprintf(f, '\tDimensions ntax=%d;\n', length(leafnames_nr));
fprintf(f, '\tTaxlabels\n');
for i = 1 : length(leafnames_nr)
    fprintf(f, '\t\t%s\n', name{leafnames_nr(i)});
end
fprintf(f, ';\n');
fprintf(f, 'End;\n');
fprintf(f, 'Begin trees;\n');
fprintf(f, '\tTranslate\n');
for i = 1 : length(leafnames_nr)
    print_str = sprintf('\t\t%d %s', leafnames_nr(i),name{leafnames_nr(i)});
    if i < length(leafnames_nr)
        print_str = [print_str ','];
    end    
    fprintf(f, '\t\t%s\n', print_str);
end
fprintf(f, ';\n');



fprintf(f, 'tree TREE1 = [&R] %s\n', regexprep(print_tree_sanity, ')n(\d*);',');'));
fprintf(f, 'End;');
fclose(f);

