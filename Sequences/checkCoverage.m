% check statistics for the per site coverage for each isolate
clear
% get all log files in the vcf folder
vcf_folders = dir('vcf/*');
files = cell(0,0);
batch = cell(0,0);
for i = 1 : length(vcf_folders)
    if vcf_folders(i).isdir && length(vcf_folders(i).name)>3
        % get all vcf files
        vcf_files_in_folder = dir(['vcf/' vcf_folders(i).name '/*.log']);
        for j = 1 : length(vcf_files_in_folder)
            files{end+1,1} = ['vcf/' vcf_folders(i).name...
                '/' vcf_files_in_folder(j).name];
            batch{end+1,1} = vcf_folders(i).name;
        end
    end    
end

% define the segments
segments = { 'HA'    'M'    'NA'    'NP'    'NS'    'PA'    'PB1'    'PB2'};

% get the minimum, median and the 10% quantile for per site coverage for
% each segment individually
for i = 1 : length(files)
    f = fopen(files{i});
    t = textscan(f, '%s\t%d\t%d'); 
    disp(files{i})
    curr_pos = 0;
    for j = 1 : length(segments)
        indices = find(ismember(t{1,1},segments{j}));
        coverage = t{1,3}(indices);
        if isempty(indices)
            stats(i,j) = 0;
        else
            stats(i,j) = quantile(coverage,0.2);
        end
            
    end
    fclose('all');
end

%%
% compare files with the same id (which run had higher coverage)
for i = 1 : length(files)
    tmp = strsplit(files{i},'-');tmp = strsplit(tmp{1},'/');
	tmp = strsplit(tmp{end},'_');
	tmp = strsplit(tmp{1},'.');
    id{i} = tmp{1};
end

% get the unique ID's
uni_id = unique(id);
is_NO = false(length(uni_id),1);
for i = 1 : length(uni_id)
    same_id = find(ismember(id,uni_id{i}));    
    fname_indices{i} = same_id;
    for j = 1 : length(same_id)
        if ~isempty(strfind(files{same_id(j)},'N'))
            is_NO(i) = true;
        end
        statistic{i,1}(j,:) = stats(same_id(j),:);
    end
    % get the minimal 20% quantile coverage of any segment 
    min_cov = min(statistic{i,1}');
    [max_cov,max_ind] = max(min_cov);
    max_min_cov(i) = max_cov; 
    mean_min_cov(i) = mean(statistic{i,1}(max_ind,:)); 
    
    % check that if the the "chosen" run has at least half of the segments
    if sum(statistic{i,1}(max_ind,1)>100)==1 && sum(statistic{i,1}(max_ind,1:end)>100)>2
        has_four(i) = true;
    else
        has_four(i) = false;
    end
        
        
    use_files{i} = files{same_id(max_ind)};
end


%%

% write all the files that should be used to a table
f = fopen('use_sequences.tsv', 'w');
use_files_ind = find(has_four);
for i = 1 : length(use_files_ind)
    fprintf(f, '%s\n', use_files{use_files_ind(i)});
end
fclose(f);

% get the files that will not be used, order them first by if they miss any
% of the first 3 segments. 

% get all the isolates that "miss" any of the segment
dont_use_files = find(max_min_cov<100);
for i = 1 : length(dont_use_files)
    if sum(statistic{dont_use_files(i)}(1:3)<100)>0
        misses_three(i) = true;
    else
        misses_three(i) = false;
    end
end

% truncate the vectors
t_mean_min_cov = mean_min_cov(dont_use_files);
t_max_min_cov = max_min_cov(dont_use_files);

% sort the 
[~,s] = sort(t_mean_min_cov);
dont_use_files(s)

f = fopen('dont_use_sequences.tsv', 'w');
g = fopen('re_sequence.tsv', 'w');
use_files_ind = find(t_max_min_cov>=100);
fprintf(g, 'isolate\tmean 0.2 q cov\n');
fprintf(f, 'isolate\tmean 0.2 q cov\tmin 0.2 q cov\tmisses one of three\tseq runs\n');
for i = length(dont_use_files) : -1 : 1
    tmp = find(ismember(id, uni_id{dont_use_files(s(i))}));
    fprintf(f, '%s\t%f\t%f\t%s\t%d\n', uni_id{dont_use_files(s(i))},... 
    t_mean_min_cov(s(i)),t_max_min_cov(s(i)), num2str(misses_three(s(i))),...
    length(tmp));
    if length(tmp)==1 && misses_three(s(i))
        fprintf(g, '%s\t%f\n', uni_id{dont_use_files(s(i))},... 
            t_mean_min_cov(s(i)));       
    end
end
fclose(f);fclose(g);





%% get the indices of the NO's (different methods)
close all
is_NO_ind = find(is_NO);
l_type = {'-',':','.'};
for i = 1 : length(is_NO_ind)
    NO_stats = statistic{is_NO_ind(i)};
    fnames = files(fname_indices{is_NO_ind(i)});
    m_cov = min(statistic{is_NO_ind(i),1}');
    [~,max_ind] = max(m_cov);
    disp(fnames(max_ind))    
    
    subplot(6,3,i)
    for j = 1 : length(fnames)
        do_plot = false;
        p_s = 0;        
        if ~isempty(strfind(fnames{j},'N_'))
            col = 'r';
            do_plot = true;
        elseif ~isempty(strfind(fnames{j},'NO_'))
            col = 'y';
            do_plot = true;
        elseif ~isempty(strfind(fnames{j},'O_'))
            col = 'g';
            f = fopen(fnames{j}); t = textscan(f, '%s\t%d\t%d'); 
            do_plot = true;
        end
        if do_plot        
            f = fopen(fnames{j}); t = textscan(f, '%s\t%d\t%d');                 
            p_s = p_s+1;
            plot(t{1,3}(1:1:end), col, 'LineStyle', l_type{p_s}); fclose(f);
        end        
        
        hold on
    end
    title(id{is_NO_ind(i)})
end
print('Coverage','-dpdf')



mins = min(stats');
[a,b] = sort(mins);
fb = files(b);

