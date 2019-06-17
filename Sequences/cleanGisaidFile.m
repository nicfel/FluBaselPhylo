% change the file from gisaid, such that there are no unneeded linebreaks
clear
f = fopen('gisaid/gisaid_full.fasta');
g = fopen('gisaid/gisaid_cleaned.fasta', 'w');
% print the first line
line = regexprep(fgets(f), '\n','');
fprintf(g, '%s', line);
while ~feof(f)
    line = regexprep(fgets(f), '\n','');
    if strcmp(line(1), '/');
        fprintf(g, '%s', line);
    else
        fprintf(g, '\n%s', line);
    end   
end
fclose('all');

