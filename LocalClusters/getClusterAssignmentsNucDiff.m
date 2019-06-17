% calculate which basel sequences are in the same cluster using parsimony
% ancestral state reconstruction
clear; fclose('all');


% read in all sampling times
seqs_file = fastaread('../Clusters/combined/HA.fasta');
bas_name = cell(0,0);
for j = 1 : length(seqs_file)
    if ~isempty(strfind(seqs_file(j).Header,'Basel'));
        tmp = strsplit(seqs_file(j).Header,'|');
        tmp = strsplit(tmp{1}, '/');
        bas_name{j} = tmp{end};
    else
        bas_name{j} = 'lala';
    end
end

% get all tree files in the combined folder
tree_files = dir('constcoalnucdiff/combined/*.trees');

% get all the tip Locations
f = fopen('tipLocations.csv'); fgets(f); c = 1;
while ~feof(f)
    line = strsplit(fgets(f),','); 
    id{c,1} = line{1};
    if strcmp(line{5},'basel')
        location{c,1} = strtrim(line{5});c = c + 1;    
    else
        location{c,1} = strtrim(line{6});c = c + 1;    
    end
end
fclose(f);

% get the unique locations
uni_loc = {'basel' 'notbasel'};

% add and unknown location
uni_loc{end+1} = 'unknown';

isBasel = find(ismember(uni_loc, 'basel'));
isUnkown = find(ismember(uni_loc, 'unknown'));

system('rm -r typednucdiff');
system('mkdir typednucdiff');
system('mkdir typednucdiff/combined');

for i = 1 : length(tree_files)
    disp(i)
    % read in the trees with the tip labels
    f = fopen(['constcoalnucdiff/combined/' tree_files(i).name]); clear name loc
    g = fopen(['typednucdiff/combined/' tree_files(i).name], 'w'); 
    trees = cell(0,0);
    tree_name = cell(0,0);
    phytrees = cell(0,0);
    print_line = true;
    
    name = cell(0,0);
    loc = zeros(0,0);
    altname = cell(0,0);
    while ~feof(f)
        full_line = fgets(f);
        line = strsplit(strtrim(full_line));
        if length(line)==2
            val = str2double(line{1});
            if ~isnan(val)
                tmp = strrep(line{2},'''','');
                name{val,1} = strrep(tmp,',','');
%                 ind = find(ismember(id, name{val}));
%                 loc(val,1) = find(ismember(uni_loc, location{ind}));
                altname{val,1} = name{val,1};
                if ~isempty(strfind(name{val,1}, 'Basel'))
                    loc(val,1) = (1);
                    tmp = strsplit(name{val}, '/');
                    tmp = strsplit(tmp{end}, '|');
                    altname{val,1} = tmp{1};
                    
                    ind = find(ismember(bas_name,tmp{1}));
                    tmp = strsplit(seqs_file(ind).Header, '|');
                    tmp2 = strsplit(strtrim(tmp{4}),'-');
                    deztime = (datenum(tmp{4},'yyyy-mm-dd')- datenum(tmp2{1},'yyyy'))...
                        /(datenum(num2str(str2double(tmp2{1})+1),'yyyy')-datenum(tmp2{1},'yyyy'))...
                        +str2double(tmp2{1});
                    samp_time(val) = deztime;                    
                else
                    samp_time(val) = 0; 
                    loc(val,1) = 2;
                end
                    
            end
            fprintf(g,'%s',full_line);
        elseif length(line)==4
            nr_nodes = strfind(line{end}, '):');
            trees{end+1,1} = sprintf(strrep(line{end},'):',')n%d:'),1:length(nr_nodes));      
            trees{end,1} = strrep(trees{end,1}, ');',sprintf(')n%d;',length(nr_nodes)+1));
            tree_name{end+1,1} = sprintf('%s %s %s',line{1}, line{2}, line{3});
            print_line = false;
        elseif print_line
            fprintf(g,'%s',full_line);
        end            
    end
    fclose(f);
    
    
    % get for each tree which sequences are in the same cluster
    [mat, cl_size, cl_source, baselNames, typetrees] = getInCluster(trees, loc, isBasel, uni_loc, isUnkown, samp_time);
    leafnames{i} = altname(baselNames);
    in_clusters{i} = mat(baselNames,baselNames);
    cluster_size{i} = cl_size(baselNames,:);
    cluster_source{i} = cl_source;
    cluster_names{i} = altname;
    
    
    % print the typetrees to file
    for j = 1 : length(typetrees)
        fprintf(g, '%s %s\n',tree_name{j}, regexprep(typetrees{j}, 'n(\d*)', ''));
    end
    fprintf(g, 'END;\n');
    fclose(g);
   
    
end


%% write the output to file
f = fopen('localclustersnucdiff/parsimony_clusters.csv', 'w');
for i = 1:length(in_clusters)
    for a = 2 : size(in_clusters{i})
        for b = a+1 : size(in_clusters{i})
            fprintf(f, '%s,%s,%f\n',leafnames{i}{a}, leafnames{i}{b}, in_clusters{i}(a,b));
        end
    end    
end
fclose(f);

f = fopen('localclustersnucdiff/parsimony_cluster_size.csv', 'w');
for i = 1: length(cluster_size)
    for a = 1 : size(cluster_size{i},1)
        fprintf(f, '%s',leafnames{i}{a});
        for b = 1 : size(cluster_size{i},2)
            fprintf(f, ',%d',cluster_size{i}(a,b));
        end
        fprintf(f, '\n');
    end    
end
fclose(f);


f = fopen('localclustersnucdiff/parsimony_cluster_source.csv', 'w');
% print the location labels in the header
line = 'rep';
for i = 1 : length(uni_loc)
    line = [line ',' uni_loc{i}];
end
fprintf(f, '%s\n', strrep(line, 'rep,', ''));
for i = 1 : length(cluster_source)
    for a = 1 : length(cluster_source{i})
        line = 'rep';
        for b = 1 : size(cluster_source{i}{a},1)
            line = [line sprintf(',%d:',cluster_source{i}{a}{b,2})];
            line = [line sprintf('{%s}',strtrim(sprintf('%d ',cluster_source{i}{a}{b,1})))];
        end
        fprintf(f, '%s\n', strrep(line, 'rep,', ''));
    end    
    fprintf(f, '-------\n');
end
fclose(f);

f = fopen('localclustersnucdiff/parsimony_cluster_membership.csv', 'w');
% print the location labels in the header
line = 'rep';
for i = 1 : length(uni_loc)
    line = [line ',' uni_loc{i}];
end
fprintf(f, '%s\n', strrep(line, 'rep,', ''));
for i = 1 : length(cluster_source)
    for a = 1 : length(cluster_source{i})
        line = 'rep';
        % iterate over clusters
        for b = 1 : size(cluster_source{i}{a},1)
            % iterate over members
            tmp_line = 'rep';
            for c = 1 : length(cluster_source{i}{a}{b,4})
                tmp_line = [tmp_line sprintf('|%s',cluster_names{i}{cluster_source{i}{a}{b,4}(c)} )];
            end
            tmp_line = strrep(tmp_line, 'rep|', '');
            line = [line ',' tmp_line sprintf(':%f:{%s}', cluster_source{i}{a}{b,3},...
                strtrim(sprintf('%d ',cluster_source{i}{a}{b,1})))];
        end
        fprintf(f, '%s\n', strrep(line, 'rep,', ''));
    end    
    fprintf(f, '-------\n');
end
fclose(f);



