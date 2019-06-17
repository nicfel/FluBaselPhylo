clear
vcf_folders = dir('vcf/*');
files = cell(0,0);
batch = cell(0,0);
for i = 1 : length(vcf_folders)
    if vcf_folders(i).isdir && length(vcf_folders(i).name)>3        
        % get all vcf files
        vcf_files_in_folder = dir(['vcf/' vcf_folders(i).name '/*.vcf']);
        for j = 1 : length(vcf_files_in_folder)
            files{end+1,1} = [vcf_files_in_folder(j).name];
            batch{end+1,1} = vcf_folders(i).name;
        end
    end    
end

uni_files = unique(files);
dublicated = cell(0,0);
for i = 1 : length(files)
    f = fopen(['vcf/' batch{i} '/' files{i}]);
    counter = 1;
    name = strsplit(files{i},'.');
    name = strsplit(name{1},'-');
    segment_name = '';

    while ~feof(f)
        line = fgets(f);
        if ~strcmp(line(1),'#')
            tmp = strsplit(line);
            if ~strcmp(segment_name,tmp{1})
                counter=1;
                segment_name = tmp{1};
            end
            af = regexp(tmp{8},'AF=(\d*).(\d*)','match');
            af = strsplit(af{1},'=');
            Data.([ batch{i} '_' name{1}]).(tmp{1}).af(counter,1) = str2double(af{2});
            Data.([ batch{i} '_' name{1}]).(tmp{1}).pos(counter,1) = str2double(tmp{2});
            Data.([ batch{i} '_' name{1}]).(tmp{1}).cal{counter,1} = tmp{5};
            Data.([ batch{i} '_' name{1}]).(tmp{1}).qual(counter,1) = str2double(tmp{6});
            counter = counter+1;
        end
    end
    fclose(f)
end
    
dublicated = unique(dublicated);
%%

for i = 1 : length(dublicated)
    name = strsplit(dublicated{i},'.');
    name = strsplit(name{1},'-');

    name1 = ['Batch1_' name{1}];
    name2 = ['Batch2_' name{1}];
    segments = fieldnames(Data.(name1));
    for j = 1:length(segments)
        if length(Data.(name1).(segments{j}).pos) ~= length(Data.(name2).(segments{j}).pos)
            fprintf('length: %s %s %d %d\n', dublicated{i}, segments{j},...
                length(Data.(name1).(segments{j}).pos),...
                length(Data.(name2).(segments{j}).pos));
            
            p1 = Data.(name1).(segments{j}).pos;
            p2 = Data.(name2).(segments{j}).pos;
            
            % get the differences in position
            p = intersect(p1,p2);
            
            i1 = find(ismember(p,p1));
            i2 = find(ismember(p,p2));
            
            p1(i1) = [];
            p2(i2) = []; 
            
            indices1 = find(ismember(Data.(name1).(segments{j}).pos,p1));
            indices2 = find(ismember(Data.(name2).(segments{j}).pos,p2));
            fprintf('\t1: %s\n',sprintf('%d ', Data.(name1).(segments{j}).af(indices1)));
            fprintf('\t2: %s\n',sprintf('%d ', Data.(name2).(segments{j}).af(indices2)));
           
%             diff = find(
        else
            if Data.(name1).(segments{j}).pos ~= Data.(name2).(segments{j}).pos
                fprintf('position: %s %s\n', dublicated{i}, segments{j})
            end
        end
    end
end
