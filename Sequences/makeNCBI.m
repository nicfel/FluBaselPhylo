clear 
% read in the metadata
f = fopen('../NonSequenceData/Master_table_processed.csv');
c = 1; 
% skip first line (headers)
fgets(f);
while ~feof(f)
    % replace empty cells with point
    line = strsplit(fgets(f),',');
    if strcmp(strtrim(line{29}), '0')
        MetaData.id{c,1} = line{2}; 
        MetaData.age(c,1) = str2double(line{12});  
        gender = strrep(line{14}, 'W', 'female');
        gender = strrep(gender, 'M', 'male');
        MetaData.gender{c,1} = gender;  
        MetaData.code{c,1} = line{1};
        datestr = strrep(line{6}, '.17','.2017');
        datestr = strrep(datestr, '.16','.2016');
        MetaData.collection{c,1} = datetime(datestr,'InputFormat','dd.MM.yy');
        c = c+1;    
    end
end
fclose(f);

f = fopen('ID_Mapping_influenza.csv');fgets(f);fgets(f);fgets(f);fgets(f);
c = 1;
while ~feof(f)
    line = strsplit(strtrim(fgets(f)), ';');
    tmp = strsplit(line{1}, '-');
    % check if id already is there
    ind = [];
    if c>1
        ind = find(ismember(Map.id, tmp{1}));
    end
    if isempty(ind)
        Map.id{c,1} = tmp{1};
        Map.newame{c,1} = strrep(line{2}, 'NMB', 'USB');
        Map.number{c,1} = strrep(line{2}, 'NMB00', '');
        c = c+1;
    end
end
fclose(f);


%%


segments = {'HA' 'M' 'NA' 'NP' 'NS' 'PA' 'PB1' 'PB2'};
segment_number = [4, 7, 6, 5, 8, 3, 2, 1];

f = fopen('NCBI_seq_map.tsv', 'w');
% print header of csv file
fprintf(f, 'SeqID\tOrganism\tcountry\tcollection-date\thost\tserotype\tPassage History\tisolation-source\tsegment\tnote\n');

c = 1;clear Data

for i = 1 : length(segments)
    % get all the Basel City HA Consensus sequences 
    clear fasta
    fasta = fastaread(['consensus/' segments{i} '.fasta']); 
    for j = 1 : length(fasta)
        id = fasta(j).Header;
        ind = find(ismember(MetaData.id, id));
        % check if the sequences only contain N's
        l = length(unique(fasta(j).Sequence));
        if ~isempty(ind) && l>1
            % remove all trailing and leadings N's for things to be
            % accepted by NCBI
            for k = length(fasta(j).Sequence) : -1 :1 
                if strcmp(fasta(j).Sequence(k), 'N')
                    fasta(j).Sequence(k) = [];
                else
                    break;
                end
            end
            while strcmp(fasta(j).Sequence(1), 'N')
                fasta(j).Sequence(1) = [];
            end

            
            mapid = find(ismember(Map.id, id));
            if isempty(mapid)
                error(sprintf('%s',id));
            end
            SeqID = sprintf('%s.%d', Map.newame{mapid}, segment_number(i));
            % get gender and age
            gender = sprintf('; %s; %d years', MetaData.gender{ind}, floor(MetaData.age(ind)));
            gender = strrep(gender, '; .','');
            disp(mapid)
            
            year = strsplit(char(MetaData.collection{ind}), '-');
            
            
            fprintf(f, '%s\tInfluenza A virus (A/Switzerland/%s/%s(H3N2))\tSwitzerland: Basel\t%s\tHomo sapiens%s\tH3N2\tDirect\tNasopharyngeal\t%d\t\n',...
                SeqID, Map.number{mapid},year{3},... 
                MetaData.collection{ind}, gender, segment_number(i));
            Data(c) = fasta(j);
            Data(c).Header = SeqID;
            c = c+1;
        end
    end
end

fclose(f);

delete('ncbi.fasta');
fastawrite( 'ncbi.fasta',Data)