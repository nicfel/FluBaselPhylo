% builds a dataset with GISAID and Basel Sequences
clear
filename = 'gisaid/gisaid_cleaned.fasta';
% read the GISAID background sequences
fasta = fastaread(filename);

% get rid unneeded line breaks that cause weired sequences
for i = 1 : length(fasta)
    if ~isempty(strfind(fasta(i).Sequence, '|'))
        tmp = strfind(fasta(i).Sequence, '|');
        fasta(i).Header = [fasta(i).Header fasta(i).Sequence(1:tmp(end))];
        fasta(i).Sequence = fasta(i).Sequence(tmp(end)+1:end);
    end
end

% get the mapping of locations
f = fopen('../NonSequenceData/geo_synonyms.tsv');
c = 1;
fgets(f);
while ~feof(f)
    line = strsplit(strtrim(fgets(f)));
    map.label{c} = lower(line{1});
    map.country{c} = lower(line{2});
    map.location{c} = lower(line{3});
    c = c+1;
end
fclose(f);

% get the mapping of regions
f = fopen('../NonSequenceData/geo_regions.tsv');
c = 1;
fgets(f);
while ~feof(f)
    line = strsplit(strtrim(fgets(f)));
    regions.country{c} = lower(line{1});
    regions.regions{c} = lower(line{2});
    regions.regions2{c} = lower(line{3});
    c = c+1;
end
fclose(f);

% list that keeps track of locations that were not assigned correctly
unknown = cell(0,0);
GISAID_accession = cell(0,0);

segs = {'HA'    'MP'    'NA'    'NP'    'NS'    'PA'    'PB1'    'PB2'};
c = ones(length(segs),1);
for i = 1 : length(fasta)
    tmp = strsplit(fasta(i).Header,'|');
    if length(tmp)>12
        % get the segment of the sequence
        segment = strtrim(tmp{13});
        % get the regional label
        loc_init = strsplit(lower(tmp{1}),'/');
        loc_init = strrep(loc_init{2}, ' ','');
        loc_init = strsplit(loc_init, '-');
        loc_init = loc_init{1};
        % look for the index of the regional label
        ind = find(ismember(map.label,loc_init));
        if ~isempty(ind)
            location = map.country{ind};
            ind2 = find(ismember(regions.country,location));          
            if ~isempty(ind2) 
                % get the GISAID accession number and the segment
                GISAID_accession{end+1,1} = [strtrim(tmp{2}) ' ' segment];                
                
                % build
                new_header = [strrep(tmp{1}, ' ', '') '|flu|' strrep(strtrim(tmp{2}) ,'_ISL_','')];
                new_header = [new_header '|' strtrim(tmp{6}) '|' regions.regions{ind2}];
                new_header = [new_header '|' regions.country{ind2} '|' location];
                new_header = [new_header '|' regions.regions2{ind2} '|?'];
                lab=strrep(strtrim(tmp{11}), ' ','_');
                if isempty(lab)
                    new_header = [new_header '|?|?|?'];
                else
                    new_header = [new_header '|' lab '|?|?'];
                end
                dt = strsplit(strtrim(tmp{6}),'-');
                seg_ind = find(ismember(segs,segment));
                Data.(segment)(c(seg_ind)) = fasta(i);
                Data.(segment)(c(seg_ind)).Header = new_header;
                c(seg_ind) = c(seg_ind)+1;
            else
                disp(location)
            end        
        else
            unknown{end+1,1} = loc_init;
        end
    else
        error('inclomplete seqeunce information');
    end
end

%% remove all dublicates
for i = 1 : length(segs)
    name = cell(length(Data.(segs{i})),1);
    for j = 1 : length(Data.(segs{i}))
        name{j} = Data.(segs{i})(j).Header;
    end
    [~,a,~]=unique(name);
    Data.(segs{i}) = Data.(segs{i})(a);
end
%%
for i = 1 : length(segs)
    delete(['gisaid/gisaid_' segs{i} '.fasta']);
    fastawrite(['gisaid/gisaid_' segs{i} '.fasta'],Data.(segs{i}));
end

%% 

% read in the sequencing metadata
f = fopen('PCR_times.csv');
i = 1; clear index time
% skip the first line
fgets(f);
while ~feof(f)
    tmp = strsplit(fgets(f),',');
    index{i,1} = tmp{1};    
    % get the time
    tmp2 = strsplit(tmp{2},'.');
    if length(tmp2) == 3
        year = strtrim(['20' tmp2{3}]);
        newtime = [year '-' tmp2{2} '-' tmp2{1}];
    else
        newtime = tmp2;
    end    
    time{i,1} = newtime;
    i = i+1;
end

system('rm -r tmp');
system('mkdir tmp');

segments = {'HA' 'M' 'NA' 'NP' 'NS' 'PA' 'PB1' 'PB2'};
for i = 1 : length(segments)
    % get all the Basel City HA Consensus sequences 
    clear fasta
    fasta = fastaread(['consensus/' segments{i} '.fasta']);
    
    c = 1; clear Data
    
    % build the names analoge to gisaid names with regions etc. to later
    % run it through the augur pipeline
    for j = 1 : length(fasta)

        tmp = strsplit(fasta(j).Header,'-');
        id = tmp{1};clear tmp
        ind = find(ismember(index,id));
        if length(ind)>0
            if length(ind)>1
                fprintf('double %s\n', id)
            end

            % build the new header with metadata
            newHeader = strcat('A/CH/Basel/', id, '|flu|', id, '|', time{ind(1)}, '|europe|switzerland|basel|this_study|cell|uni_hospital_basel|00|?');
            Data(c) = fasta(j);
            Data(c).Header = newHeader;    
            c = c+1;
        else
            fprintf('none %s\n', id)
        end
    end
    delete(['tmp/' segments{i} '.fasta']);
    fastawrite(['tmp/' segments{i} '.fasta'], Data);
end


%% make the augur input files
HA_bas = fastaread('tmp/HA.fasta');
HA_gis = fastaread('gisaid/gisaid_HA.fasta');
HA_anc = fastaread('gisaid/ancestral_HA.fasta');


delete('../TimeTree/augur/builds/flu/h3n2_ha.fasta');
fastawrite('../TimeTree/augur/builds/flu/h3n2_ha.fasta', HA_bas);
fastawrite('../TimeTree/augur/builds/flu/h3n2_ha.fasta', HA_gis);
fastawrite('../TimeTree/augur/builds/flu/h3n2_ha.fasta', HA_anc);
