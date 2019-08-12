% combine euler age mixing output
clear
% check which interval sizes and cutoff were used
all_files = dir('groupmixing/group_mixing_pairs*.csv');
int_size = zeros(0,0);
cut_offs = zeros(0,0);
for i  = 1 : length(all_files)
    tmp = strsplit(all_files(i).name, '_');
    cut_offs(i) = str2double(tmp{4});
    int_size(i) = str2double(tmp{5});
end
% get the unique ones
int_size = unique(int_size);
cut_offs = unique(cut_offs);

g = fopen('out/file_values_group_mixing.csv', 'w');
fprintf(g, 'filename,lower,upper,interval_size\n');

% file count
f_count = 1;

p_correction = ((6^2-6)/2+6)*2;


% for a = 1 : length(cut_offs)
%     for b = 1 : length(int_size)
%         clear vals
%         input1 = dir(sprintf('groupmixing/group_mixing_%d_%d_*.csv', cut_offs(a), int_size(b)));
% 
%         clear linestart
%         % loop over all files
%         for i = 1 : length(input1)
%             disp(input1(i).name)
%             f = fopen(['groupmixing/' input1(i).name]);        
%             t = textscan(f, '%s %s %f %f','Delimiter',',','HeaderLines', 3 );
%             fclose(f);
% 
%             
%             if i==1
%                 for j = 1 : length(t{1})
%                     linestart{j} = [t{1}{j} ',' t{2}{j}];
%                 end
%                 vals(:,1) = t{3};
%                 vals(:,2) = t{4};
%             else
%                 vals(:,1) = vals(:,1) + t{3};
%                 vals(:,2) = vals(:,2) + t{4};
%             end
% 
% 
% 
%         end
% 
%         % read in the cutoff values
%         f = fopen(['groupmixing/' input1(1).name]); 
%         fprintf(g, '%s',sprintf('out/group_mixing_%d.csv',f_count));
%         tmp = strsplit(strtrim(fgets(f)));
%         fprintf(g, ',%s',tmp{end});
%         tmp = strsplit(strtrim(fgets(f)));
%         fprintf(g, ',%s',tmp{end});
%         fprintf(g, ',%d\n',int_size(b));
%         fclose(f);
% 
% 
% 
%     %%
%         old_vals = vals;
%         clear prob perc
% 
%         % take the average
%         for k = 1 : size(vals,1)
%             prob(k,1) = 1 - vals(k,1)/(length(input1)*100000);
%             prob(k,2) = 1 - vals(k,2)/(length(input1)*100000);
%         end
% 
% 
%         % scale the percentile to log scale
%         for k = 1 : size(prob,1)
%             if prob(k,1)>prob(k,2)
%                 perc(k,1) = log10(min(1,prob(k,2)*p_correction));
%             else
%                 perc(k,1) = -log10(min(1,prob(k,1)*p_correction));
%             end
%         end   
% 
% 
%         % print to file again
%         f = fopen(sprintf('out/group_mixing_%d.csv',f_count),'w');
%         % print header
%         fprintf(f,'from,to,percentile\n');
%         for k = 1 : size(prob,1)
%             tmp = strsplit(linestart{k}, ',');
%             fprintf(f, '%s,%s,%s\n',tmp{1}, tmp{2}, num2str(perc(k)));
%             if ~strcmp(tmp{1}, tmp{2})
%                 fprintf(f, '%s,%s,%s\n',tmp{2}, tmp{1}, num2str(perc(k)));
%             end
%         end
%         f_count = f_count + 1;
%         fclose(f);   
%     end
% end
% fclose('all');

%% do the same based on pairs instead of unique patients
% file count
f_count = 1;
g = fopen('out/file_values_group_mixing_pairs.csv', 'w');
fprintf(g, 'filename,lower,upper,interval_size\n');


for a = 1 : length(cut_offs)
    for b = 1 : length(int_size)
        clear vals
        input1 = dir(sprintf('groupmixing/group_mixing_pairs_%d_%d_*.csv', cut_offs(a), int_size(b)));

        clear linestart
        % loop over all files
        for i = 1 : length(input1)
            disp(input1(i).name)
            f = fopen(['groupmixing/' input1(i).name]);        
            t = textscan(f, '%s %s %f %f','Delimiter',',','HeaderLines', 3 );
            fclose(f);            
            if i==1
                for j = 1 : length(t{1})
                    linestart{j} = [t{1}{j} ',' t{2}{j}];
                end
                vals(:,1) = t{3};
                vals(:,2) = t{4};
            else
                vals(:,1) = vals(:,1) + t{3};
                vals(:,2) = vals(:,2) + t{4};
            end
        end

        % read in the cutoff values
        f = fopen(['groupmixing/' input1(1).name]); 
        fprintf(g, '%s',sprintf('out/group_mixing_pairs_%d.csv',f_count));
        tmp = strsplit(strtrim(fgets(f)));
        fprintf(g, ',%s',tmp{end});
        tmp = strsplit(strtrim(fgets(f)));
        fprintf(g, ',%s',tmp{end});
        fprintf(g, ',%d\n',int_size(b));
        fclose(f);



    %%
        old_vals = vals;
        clear prob perc

        % take the average
        for k = 1 : size(vals,1)
            prob(k,1) = 1 - vals(k,1)/(length(input1)*100000);
            prob(k,2) = 1 - vals(k,2)/(length(input1)*100000);
        end


        % scale the percentile to log scale
        for k = 1 : size(prob,1)
            if prob(k,1)>prob(k,2)
                perc(k,1) = log10(min(prob(k,2)*p_correction,1));
            else
                perc(k,1) = -log10(min(prob(k,1)*p_correction,1));
            end
        end   


        % print to file again
        f = fopen(sprintf('out/group_mixing_pairs_%d.csv',f_count),'w');
        % print header
        fprintf(f,'from,to,percentile\n');
        for k = 1 : size(prob,1)
            fprintf(f, '%s,%s\n',linestart{k}, num2str(perc(k)));
        end
        f_count = f_count + 1;
        fclose(f);   
    end
end
fclose('all');
