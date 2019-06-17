function [in_cluster, cluster_size, cluster_source, basel_names, typetrees] = getInCluster(trees, loc, isBasel, loc_names, isUnknown, samp_time)
% get the unique locations
uni_loc = unique(loc);

% get which basel sequences are grouped together
isInBasel = find(loc==isBasel);

in_cluster_mat = zeros(length(loc), length(loc));
cluster_size = zeros(length(loc), 1);


clear isInBasel

wait = waitbar(0,'Please wait...');

% for each tree, get the cluster assignment of each leaf
for i = 1 : length(trees)
    waitbar(i/length(trees))

    % read in the tree as a phytree object
    ptree = phytreeread(trees{i});
    nodenames = get(ptree, 'nodenames');
    
    
    % get the connectivity matrix and the distances between nodes
    mat = getmatrix(ptree);
    dist = pdist(ptree, 'Squareform',true,'Nodes','all');

    
    % initialize the location vector
    location = cell(length(nodenames),1);
    visited = false(length(nodenames),1);
    dist_from_basel_sample = zeros(length(nodenames),1);
    node_date = zeros(length(nodenames),1);

    for j = 1 : length(loc)
        ind = str2double(nodenames{j});
        location{j} = loc(ind);        
        visited(j) = true;
        if loc(ind) == isBasel
            dist_from_basel_sample(j) = 0;
            node_date(j) = samp_time(ind);
        else
            dist_from_basel_sample(j) = Inf;
            node_date(j) = Inf;
        end
    end
    
    % upwards path of the parsimony calling    
    not_visited = find(~visited);
    while ~isempty(not_visited)
        not_visited = find(~visited);
        for j = 1 : length(not_visited);
            children = find(mat(not_visited(j),:));
            if ~isempty(location{children(1)}) && ~isempty(location{children(2)})
                % get the distance of the new sample from a basel sequence
                pdist_sample = dist(not_visited(j), children);
                dist_from_basel_sample(not_visited(j)) = min(dist_from_basel_sample(children) + pdist_sample');
                node_date(not_visited(j)) = min(node_date(children) - pdist_sample');

                int = intersect(location{children(1)}, location{children(2)});
                if isempty(int)
                    location{not_visited(j)} = [location{children(1)}, location{children(2)}];
                else
                    location{not_visited(j)} = int;
                end
                visited(not_visited(j)) = true; 
                
%                 ind_basel = find(location{not_visited(j)}==isBasel);
%                 if ~isempty(ind_basel) && length(location{not_visited(j)})>1 
%                     if dist_from_basel_sample(not_visited(j)) > 0.1
%                         location{not_visited(j)}(ind_basel) = [];
%                     end
%                 end

                
                % make sure that if the distance of the node is further than 0.1 years from
                % a basel sequence that it's location is considered unknown
%                 if node_date(not_visited(j))<2016.748634 || dist_from_basel_sample(not_visited(j)) > 0.1
                if dist_from_basel_sample(not_visited(j)) > 0.1
                    ind_basel = find(location{not_visited(j)}==isBasel);
                    if ~isempty(ind_basel)
                        % if there is more than one option, and unknown is but one of
                        % them, don't choose unknown
%                         if length(location{not_visited(j)})>1
%                             location{not_visited(j)}(ind_basel) = [];
%                         else
                            if isempty(find(ismember(location{not_visited(j)},isUnknown)))
                                location{not_visited(j)} = [location{not_visited(j)} isUnknown];
                            end
%                         end                            
                    end                        
                end                
               
%                 if length(location{not_visited(j)})>1
%                     location{not_visited(j)}(location{not_visited(j)}==isUnknown) = [];
%                 end                
            end
        end
    end
    
    % downwards calling 
    visited(length(loc)+1:end) = false;
    visited(end) = true;
    
    not_visited = find(~visited);
    while ~isempty(not_visited)
        not_visited = find(~visited);
        for j = length(not_visited) : -1 : 1
            parent = find(mat(:,not_visited(j)));            
            if sum(ismember(not_visited, parent))==0
                
                int = intersect(location{not_visited(j)}, location{parent});
                if ~isempty(int)
                    location{not_visited(j)} = int;
                end
                visited(not_visited(j)) = true;               

            end
        end
    end
    
    % choose for each node a random node location
    random_node_location = zeros(length(location),1);
    for j = 1 : length(location)
        if length(location{j})>1
            % never choose basel to get a conservative guess of the
            % ancestral location as well as for which lineages were in Basel
            tmp = location{j};
            tmp(tmp==isBasel) = [];
            if length(tmp)>1
                random_node_location(j) = randsample(tmp,1);
            else
                random_node_location(j) = tmp;
            end
        else
            random_node_location(j) = location{j}(1);
        end            
    end
    
    
    % get all node that are only in basel
    onlyBasel = false(length(visited),1);   
    for j = 1 : length(onlyBasel)
        if length(location{j})==1 && location{j}==isBasel
            onlyBasel(j) = true;
        end        
    end
    
    isInBasel = find(onlyBasel(1:length(loc)));
    
    basel_names = str2double(nodenames(isInBasel));    
    
    clustered_node = false(size(nodenames));
    clustered_node(isInBasel) = true;

    
    % for each isInBasel leaf, get all parent leaves in Basel
    baselParent = cell(length(isInBasel),1);
    for a = 1 : length(isInBasel)
        parent = find(mat(:,isInBasel(a)));
        while onlyBasel(parent)
            clustered_node(parent) = true;
            baselParent{a} = [baselParent{a} parent];
            parent = find(mat(:,parent));
        end
    end
    
    new_mat = zeros(size(in_cluster_mat));
    % get for each pair of leaves if they are in the same cluster
    for a = 1 : length(isInBasel)
        if ~isempty(baselParent{a})
            for b = a+1 : length(isInBasel)
                if ~isempty(baselParent{b})
                    int = intersect(baselParent{a}, baselParent{b});
                    if ~isempty(int)
                        new_mat(basel_names(a),basel_names(b)) = 1;
                    end
                end
            end     
        end
    end   
    in_cluster_mat = in_cluster_mat + new_mat;    
    
    % check which sequences are in the same cluster
    nr_cl = 1;
    clustered = zeros(0,0); clear cluster
    for a = 1 : length(isInBasel)
        if isempty(find(ismember(clustered,basel_names(a))))
            same_cluster = find(new_mat(basel_names(a),:));
            cluster{nr_cl} = [basel_names(a) same_cluster];
            clustered = [clustered cluster{nr_cl}];
            nr_cl = nr_cl+1;
        end
    end
    
    % get the cluster size that each leave is in
    for j = 1 : length(cluster)
        for k = 1 : length(cluster{j})
            cluster_size(cluster{j}(k),i) = length(cluster{j});
        end
    end
    
    cluster_name = zeros(size(nodenames));
    lc_name=1;
    % get the cluster size that each leave is in
    for j = 1 : length(cluster)
        for k = 1 : length(cluster{j})
            cluster_name(find(ismember(nodenames, num2str(cluster{j}(k))))) = lc_name;    
        end
        lc_name = lc_name+1;
    end
    
    
    % get the source location of each cluster
    clear source_loc
    for j = 1 : length(cluster)
        % get the first non basel parent of the basel sequence
        ind = find(ismember(basel_names, cluster{j}(1)));
        if length(baselParent{ind})>0        
            last_basel_parent = baselParent{ind}(end);
        else
            last_basel_parent = isInBasel(ind);
        end         
        
        next_parent_ind = find(mat(:,last_basel_parent));
        clear cl_leave_ind
        % get the indices of the basel leaves
        for k = 1 : length(cluster{j})
            cl_leave_ind(k) = find(ismember(nodenames,num2str(cluster{j}(k))));
        end
        
        if ~isempty(next_parent_ind)
            source_loc{j,1} = location{next_parent_ind};
            source_loc{j,2} = length(cluster{j});
            % also get the distance of that node to the most recent sample
            % from basel in that cluster
            source_loc{j,3} = max(dist(next_parent_ind,cl_leave_ind)); 
            source_loc{j,4} = sort(cluster{j});
        else
            % Should this be unknwown source instead?
            source_loc{j,1} = location{last_basel_parent};
            source_loc{j,2} = length(cluster{j}); 
            % if for some reason the source is in basel, use nan
            source_loc{j,3} = NaN;
            source_loc{j,4} = sort(cluster{j});
        end
    end
    
    cluster_source{i,1} = source_loc;  
    
    % make a new trees variable
    typetrees{i} = trees{i};
    % relable the trees string to have the locations in there
    for j = 1 : (length(nodenames)-1)/2 +1
        typetrees{i} = regexprep(typetrees{i}, ['(' nodenames{j} ':'],...
            ['(' nodenames{j} '[&loc="' loc_names{random_node_location(j)}...
            '",clust=', num2str(clustered_node(j)) ',lcn=', num2str(cluster_name(j)) ']:']);
        typetrees{i} = regexprep(typetrees{i}, [',' nodenames{j} ':'],...
            [',' nodenames{j} '[&loc="' loc_names{random_node_location(j)}...
            '",clust=', num2str(clustered_node(j)) ',lcn=', num2str(cluster_name(j)) ']:']);
    end
    
    for j = (length(nodenames)-1)/2 + 2 : length(nodenames)-1
        typetrees{i} = regexprep(typetrees{i}, [')' nodenames{j} ':'],...
            [')' nodenames{j} '[&loc="' loc_names{random_node_location(j)}...
            '",clust=', num2str(clustered_node(j)) ',lcn=', num2str(cluster_name(j)) ']:']);
    end
    % add the root location
    typetrees{i} = regexprep(typetrees{i}, [')' nodenames{end} ';'],...
        [')' nodenames{end} '[&loc="' loc_names{random_node_location(end)}...
        '",clust=', num2str(clustered_node(end)) ',lcn=', num2str(cluster_name(j)) '];']);

    
end

in_cluster = in_cluster_mat./length(trees);

close(wait)


end