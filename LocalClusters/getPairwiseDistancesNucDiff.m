% get the distances between the any two patients in the same clusters
clear
% get all tree files in the combined folder
tree_files = dir('constcoalnucdiff/combined/*.trees');

% get all the tip Locations
f = fopen('tipLocations.csv'); fgets(f); c = 1;
while ~feof(f)
    line = strsplit(fgets(f),','); 
    tmp = strsplit(line{1}, '.');
    id{c,1} = tmp{1};
    if strcmp(line{5},'basel')
        location{c,1} = strtrim(line{5});c = c + 1;    
    else
        location{c,1} = strtrim(line{6});c = c + 1;    
    end
end
fclose(f);

% get the unique locations
uni_loc = {'basel' 'notbasel'};

isBasel = find(ismember(uni_loc, 'basel'));


for i = 1 : length(tree_files)
    disp(i)
    % read in the trees with the tip labels
    f = fopen(['constcoalnucdiff/combined/' tree_files(i).name]); clear name loc
    trees = cell(0,0);
    phytrees = cell(0,0);
    while ~feof(f)
        line = strsplit(strtrim(fgets(f)));
        if length(line)==2
            val = str2double(line{1});
            if ~isnan(val)
                tmp = strrep(line{2},'''','');
                name{val,1} = strrep(tmp,',','');
                if ~isempty(strfind(name{val,1}, 'Basel'))
                    loc(val,1) = 1;
                else
                    loc(val,1) = 2;
                end
%                 loc_name = strsplit(name{val}, '|');
%                 
%                 ind = find(ismember(id, loc_name{1}));
%                 loc(val,1) = find(ismember(uni_loc, location{ind}));
            end
        end
        
        if length(line)==4
            nr_nodes = strfind(line{end}, '):');
            trees{end+1,1} = sprintf(strrep(line{end},'):',')n%d:'),1:length(nr_nodes));      
            trees{end,1} = strrep(trees{end,1}, ');',sprintf(')n%d;',length(nr_nodes)+1));
        end       
    end
    fclose(f);
    
    
    % get for each tree which sequences are in the same cluster
    [distances] = getDistances(trees, loc, isBasel);
    inBasel = find(loc==isBasel);
    leafnames{i} = cell(0,0);
    for j = 1 : length(inBasel)
        leafnames{i}{j} = name{inBasel(j)};
    end
    pairwise_distances{i} = distances(inBasel,inBasel,:);
end

%% write the output to file
f = fopen('localclustersnucdiff/pairwise_distances.csv', 'w');

% print the location labels in the header
for i = 1 : length(pairwise_distances)
    for a = 1 : size(pairwise_distances{i},1)
        for b = a+1 : size(pairwise_distances{i},2)
            tmp1 = strsplit(leafnames{i}{a},'|');
            tmp1 = strsplit(tmp1{1},'/');
            tmp1 = strsplit(tmp1{end},'.');
            
            tmp2 = strsplit(leafnames{i}{b},'|');
            tmp2 = strsplit(tmp2{1},'/');
            tmp2 = strsplit(tmp2{end},'.');
            
            line = sprintf('%s,%s',tmp1{1},tmp2{1});
            for c = 1 : size(pairwise_distances{i},3)
                line = [line sprintf(',%f',pairwise_distances{i}(a,b,c))];
            end
            fprintf(f, '%s\n', line);
        end
    end    
end
fclose(f);


%% 
all = zeros(0,0);
for i = 1 : length(pairwise_distances)
%     subplot(11,8,i)
    tmp = pairwise_distances{i}(pairwise_distances{i}>0);
    if length(tmp)>0
        all = [all; tmp];
%         ksdensity(tmp);hold on
    end
%     if length(tmp)>0
%         ksdensity(tmp);
%     end
end

ksdensity(all)

