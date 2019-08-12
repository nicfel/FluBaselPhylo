% get the number of patients per age interval
clear

% read in the metadata
f = fopen('../NonSequenceData/Master_table_processed.csv');
c = 1; 
% skip first line (headers)
fgets(f);
while ~feof(f);
    % replace empty cells with point
    line = strsplit(fgets(f),',');
    if strcmp(strtrim(line{29}), '0')
        MetaData.id{c,1} = line{2}; 
        MetaData.age(c,1) = str2double(line{12});  
        c = c+1;    
    end
end
fclose(f);



% print to file
f = fopen('out/age_distribution.csv', 'w');
% make the header
header = 'age_group,numbers';
% print header
fprintf(f,'%s\n',header);
    
upper = 10:10:100;
lower = upper - 10; 
upper(end)=100;
for a = 1 : length(upper)
    age_group_a = find(MetaData.age>=lower(a) & MetaData.age < upper(a));
    val_from = sprintf('%03d_%03d',lower(a),upper(a));
    fprintf(f, '%s,%d\n',val_from, length(age_group_a));
end
fclose(f);



