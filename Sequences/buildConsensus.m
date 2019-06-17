% build consensus sequences from the vcf and coverage files 
clear
% get the list of sequences for which consensus should be built
f = fopen('use_sequences.tsv');list = textscan(f, '%s');fclose(f);

% read in the reference sequence
% build the consensus sequences
segments = {'HA' 'M' 'NA' 'NP' 'NS' 'PA' 'PB1' 'PB2'};
system('rm -r consensus');
system('mkdir consensus');
warning('off','all')
for i = 1 : length(list{1})
    reference = fastaread('rawreads/reference.fasta');
   
    % get the filename and isoalte ID
    file_name = list{1}{i};
    isolate_id = strsplit(file_name,'-');
    isolate_id = strsplit(isolate_id{1},'/');
    isolate_id = strsplit(isolate_id{end}, '_');
    isolate_id = strsplit(isolate_id{1}, '.');
    isolate_id = isolate_id{1};
    disp(isolate_id)
    
    % read the coverage data
    f = fopen(file_name); t = textscan(f, '%s\t%d\t%d');fclose(f);
    % reaf the variant calls
    f = fopen(strrep(file_name, '.log','.vcf')); clear v
    count = 1;
    while ~feof(f)
        line = fgets(f);
        if ~strcmp(line(1), '#')
            tmp = strsplit(line);
            v{1}{count,1} = tmp{1};
            v{2}(count,1) = str2double(tmp{2});
            v{3}{count,1} = tmp{4};
            v{4}{count,1} = tmp{5};
            af = regexp(line, 'AF=(\d).(\d*)','match');
            af = strsplit(af{1}, '=');
            v{5}(count,1) = str2double(af{2});
            count = count+1;
        end
        
    end
    fclose(f);
    
    
    % do consensus calling
    for j = 1 : length(segments)
        % set N for all entries with coverage lower than 100
        indices = find(ismember(t{1,1},reference(j).Header));   
        coverage = t{1,3}(indices);
        seq = reference(j).Sequence;
        seq_cov = zeros(length(seq),1);
        seq_cov(t{1,2}(indices)) = coverage;
        low_cov = find(seq_cov<100);
        seq(low_cov) = 'N';
        
        % do variant calling
        indices_vcf = find(ismember(v{1,1},reference(j).Header));   
        af_vcf = v{1,5}(indices_vcf);
        indices_vcf(af_vcf<0.5) = [];
        positions_vcf = v{1,2}(indices_vcf);
        new_base_vcf = v{1,4}(indices_vcf);
        
%         disp(file_name);
%         disp(length(v{1,1}));

        
        for k = 1 : length(new_base_vcf)
            if ~strcmp(seq(positions_vcf(k)), 'N')
                seq(positions_vcf(k)) = new_base_vcf{k};
            end
        end
        
        new_file(1).Header = isolate_id;
        new_file(1).Sequence = seq;
        fastawrite(['consensus/' segments{j} '.fasta'], new_file);
    end
    fclose('all');
end
warning('on','all')
