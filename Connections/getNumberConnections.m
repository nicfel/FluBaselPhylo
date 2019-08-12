% get to how many individuals each individual is connected to
clear
% load pairwise distances
load('../AgeMixing/input_vals_group_mixing_nucdiff.mat')


uni_group = unique(MetaData.group);

for cut = 1 : length(c_off_lower)
    % get all value_all entries that have distances between the lower and
    % upper cutoff values
    use_value_indices = find(value_dublicated(:,4)>=c_off_lower(cut) &...
                                value_dublicated(:,4)<c_off_upper(cut));

    % only use the subset of pairs that have a distance withing the lower
    % and upper bounds
    value = value_dublicated(use_value_indices,:);

    for i = 1 : length(MetaData.id)
        ind = find(value(:,6)==i);
        for j = 1 : length(uni_group)
            nr_cons(cut,i,j) = sum(MetaData.group(value(ind,7))==j);
        end
    end
end
for cut = 1 : length(c_off_lower)
    f = fopen(sprintf('numberConnections_cuotff%d.csv', cut),'w');
    fprintf(f, 'age,group,con.group1,con.group2,con.group3,con.group4,con.group5,con.group6\n');
    for i = 1 : length(MetaData.id)
        fprintf(f, '%f,%d,%d,%d,%d,%d,%d,%d\n',MetaData.age(i), MetaData.group(i),...
            nr_cons(cut,i,1),nr_cons(cut,i,2),nr_cons(cut,i,3),...
            nr_cons(cut,i,4),nr_cons(cut,i,5),nr_cons(cut,i,6));    
    end
    fclose(f);
end

f = fopen('numberConnections.csv','w');
fprintf(f, 'age,group,cutoff1,cutoff2,cutoff3,cutoff4,cutoff5\n');
for i = 1 : length(MetaData.id)
    fprintf(f, '%f,%d,%d,%d,%d,%d,%d\n',MetaData.age(i), MetaData.group(i),...
        sum(nr_cons(1,i,:)),sum(nr_cons(2,i,:)),sum(nr_cons(3,i,:)),sum(nr_cons(4,i,:)),sum(nr_cons(5,i,:)));    
end
fclose(f);
