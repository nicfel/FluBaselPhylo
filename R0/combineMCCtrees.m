% combine the mcc trees from the bdsky analysis into one tree with
% singleton nodes for the origin
clear; fclose('all')
% get the sampling time of each sample 
fas = fastaread('../TimeTree/augur/builds/flu/h3n2_ha.fasta');
samplingtime_name = cell(0,0);
is_basel = false(0,0);
for i = 1 : length(fas)
    tmp = strsplit(fas(i).Header, '|');    
    id{i} = tmp{1};
    date = tmp{4};
    tmp2 = strsplit(date,'-');
    if length(tmp2)==3 && isempty(strfind(date, 'X'))
        deztime = (datenum(date,'yyyy-mm-dd')- datenum(tmp2{1},'yyyy'))...
                /(datenum(num2str(str2double(tmp2{1})+1),'yyyy')-datenum(tmp2{1},'yyyy'))...
                +str2double(tmp2{1});
        samplingtime(i) = deztime;
        samplingtime_name{i} = date;
        
        if ~isempty(strfind(id{i}, 'A/CH/Basel/'))
            is_basel(i) = true;
        end
    else
        samplingtime(i) = 0;
    end
end
all_names = cell(0,0);

% display the max sampling time
[~,maxind] = max(samplingtime(is_basel));
samplingtime_red = samplingtime_name(is_basel);
disp(samplingtime_red{maxind})


% get all the mcc trees
system('mkdir mcccombined')
for subset = 1 : 10
    disp(subset)
    mcc_trees = dir(['bdskynucdiff/mcc/*subs' num2str(subset) '.trees']);
    % read all the trees
    for i = 1 : length(mcc_trees)
        f = fopen(['bdskynucdiff/mcc/' mcc_trees(i).name]);
        while ~feof(f)
            full_line = fgets(f);
            line = strsplit(strtrim(full_line));
            if length(line)==2
                val = str2double(line{1});
                if ~isnan(val)
                    tmp = strrep(line{2},'''','');
                    name{val,1} = strrep(tmp,',','');
                    all_names{end+1} = name{val,1};

                    % find the sampling time of every leaf in this cl
                    tmp_id = strsplit(name{val,1},'.');
                    tmp_id = strsplit(tmp_id{1},'_');
                    true_id = tmp_id{2};
                    id_id = find(ismember(id, true_id));
                    samptime(val) = samplingtime(id_id(1));


                end
            elseif length(line)==4
                nr_nodes = strfind(line{end}, '):');
                mcc_tree{i} = line{end};
                % replace each number with the actual leafname
                for j = 1 : length(name)
                    mcc_tree{i} = strrep(mcc_tree{i}, sprintf('(%d[', j), sprintf('(%s[', name{j}) );
                    mcc_tree{i} = strrep(mcc_tree{i}, sprintf(',%d[', j), sprintf(',%s[', name{j}) );
                end

                heights_char = regexp(mcc_tree{i}, 'height=(\d*)\.(\d*),', 'match');
                clear height
                for j = 1 : length(heights_char)
                    tmp = strrep(heights_char{j}, 'height=','');
                    tmp = strrep(tmp, ',', '');
                    height(j) = str2double(tmp);
                end

                root_height(i) = max(samptime) - max(height); 
                all_heights(i) = max(samptime);
            end            
        end
    end

    %% also get all the "cluster" of size 1
    % read in the xml
    f = fopen('bdskynucdiff/xmls/bdsky_subs1_rep0.xml');
    while ~feof(f)
        line = fgets(f);
        if ~isempty(strfind(line, '<trait id="dateTrait.t:l'))
            nextline = fgets(f);
            if ~isempty(strfind(nextline, '">'))
                tmp = strsplit(strtrim(nextline), '=');
                mcc_tree{end+1} = [tmp{1} '[&]:0.0;'];
                all_names{end+1} = tmp{1};

                tmp_id = strsplit(tmp{1},'.');
                tmp_id = strsplit(tmp_id{1},'_');
                true_id = tmp_id{2};
                id_id = find(ismember(id, true_id));
                root_height(end+1) = samplingtime(id_id(1));
                all_heights(end+1) = samplingtime(id_id(1));
            end
        end
    end
    fclose(f);




    %% combine all the mcc trees including their origin
    % open the log file to get the origin


    % get the minimal root height
    min_root_heigh = min(root_height);
    origin_height = min_root_heigh-0.0001;

    sort_height = (root_height + all_heights)/2;

    % sort according to most recent sample
    [~,ij] = sort(sort_height(1:end),'descend');

    origin = 0.01;
    % combine trees into 1
    for j = 1 : length(mcc_tree)
        i = ij(j);

        tmp = mcc_tree{i};
        tmp = strrep(tmp, ']', ',known=1]');
        tmp = [strrep(tmp, ':0.0;',...
            sprintf(':%f)[&known=0]:%f', origin, abs(min_root_heigh-root_height(i))+0.0000*(length(mcc_tree)-j) ))];
        tmp = ['(' tmp];
        if j == 1
            print_tree = tmp;
        else
            print_tree = [  '(' print_tree ',' tmp ')[&known=0]:' num2str(0.0000*j) ];
        end    
    end
    print_tree = ['(' print_tree ')[&known=0]:0.1;'];
    print_tree = strrep(print_tree, '&,','&');
    f = fopen(['mcccombined/all_trees_subs' num2str(subset) '.trees'], 'w');
    fprintf(f, '#NEXUS\n');
    fprintf(f, 'Begin taxa;\n');
    fprintf(f, '\tDimensions ntax=%d;\n', length(all_names));
    fprintf(f, '\t\tTaxlabels\n');
    for i = 1 : length(all_names)
        fprintf(f, '\t\t\t%s\n', all_names{i});
    end
    fprintf(f, ';\n');
    fprintf(f, 'End;\n');
    fprintf(f, 'Begin trees;\n');
    fprintf(f,'tree TREE1 = %s\n', print_tree);
    fprintf(f, 'End;\n');
    fclose('all');
end


dasdas

% %% build the mapping of trait data
% f = fopen('../MixingPatterns/age_gender.csv');
% c = 1; 
% % skip first line (headers)
% fgets(f);
% while ~feof(f);
%     % replace empty cells with point
%     tmp_line = strrep(fgets(f),',,', ',.,');
%     line = strsplit(tmp_line,',');
%     % assign id and age label
%     MetaData.id{c} = line{1};
%     MetaData.age(c) = str2double(line{2});
%     c = c+1;    
% end
% fclose(f);

% read in the file that contains the information if people ave children in
% their household
f = fopen('../Questionnairs/HouselholdMember.csv');
c = 1; 

% skip first line (headers)
fgets(f);
while ~feof(f);
    % replace empty cells with point
    line = strsplit(fgets(f),',');

    HousMetaData.id{c} = line{1};

    age(1) = str2double(line{2});
    age(2) = str2double(line{3});

    age_ind = find(ismember(MetaData.id, HousMetaData.id{c}));

    HousMetaData.hasAnswer(c) = str2double(line{6});
    % if there is no answer and the person is below 20 to exlude
    % children mixing

    HousMetaData.isKid(c) = 0;

    if length(age_ind)>0
        HousMetaData.age(c) = MetaData.age(age_ind(1));
        if  MetaData.age(age_ind(1))<20
            HousMetaData.isKid(c) = 1;
        end
    else
        HousMetaData.age(c) = -1000;
    end

    if age(1)>0 || age(2)>0
        HousMetaData.hasKids(c) = 1;
    elseif age(1)==0 || age(2)==0
        HousMetaData.hasKids(c) = 0;
    else
        HousMetaData.hasKids(c) = -1;
    end

    if age(1)>0
        HousMetaData.hasYoungKids(c) = 1;
    elseif age(1)==0
        HousMetaData.hasYoungKids(c) = 0;
    else
        HousMetaData.hasYoungKids(c) = -1;
    end    

    if age(2)>0
        HousMetaData.hasOldKids(c) = 1;
    elseif age(2)==0
        HousMetaData.hasOldKids(c) = 0;
    else
        HousMetaData.hasOldKids(c) = -1;
    end    

    c = c+1;    
end
fclose(f);


HousMetaData.hasAnswer(HousMetaData.hasAnswer==0) = NaN;


HousMetaData.hasKids(HousMetaData.hasKids==-1) = NaN;
HousMetaData.hasYoungKids(HousMetaData.hasYoungKids==-1) = NaN;
HousMetaData.hasOldKids(HousMetaData.hasOldKids==-1) = NaN;

HousMetaData.hasKids(HousMetaData.hasKids==0) = NaN;
HousMetaData.hasYoungKids(HousMetaData.hasYoungKids==0) = NaN;
HousMetaData.hasOldKids(HousMetaData.hasOldKids==0) = NaN;



f = fopen('quest_labels.csv', 'w');
fprintf(f, 'answers,kids,young,old\n');
for i = 1 : length(all_names)
    tmp = strsplit(all_names{i}, '/');
    tmp = strsplit(tmp{end}, '.');
    ind1 = find(ismember(MetaData.id, tmp{1}));    
    ind2 = find(ismember(HousMetaData.id, tmp{1}));

    fprintf(f, '%s,%d,%d,%d,%d\n', all_names{i},...
        HousMetaData.hasAnswer(ind2(1)), HousMetaData.hasKids(ind2(1)),...
         HousMetaData.hasYoungKids(ind2(1)), HousMetaData.hasOldKids(ind2(1)));
end
fclose(f);

f = fopen('age_labels.csv', 'w');
fprintf(f, 'taxa,age\n');
for i = 1 : length(all_names)
    tmp = strsplit(all_names{i}, '/');
    tmp = strsplit(tmp{end}, '.');
    ind1 = find(ismember(MetaData.id, tmp{1}));    
    ind2 = find(ismember(HousMetaData.id, tmp{1}));

    fprintf(f, '%s,%d\n', all_names{i},...
        round(MetaData.age(ind1(1))));
end
fclose(f);
