% get the weights for the different groups
clear
% read in the number of people in each age bracket of Basel Stadt
f = fopen('../NonSequenceData/age_basel_stadt.csv');fgets(f);
while ~feof(f)
    line = strsplit(fgets(f), ',');
    tmp1 = strsplit(line{1}, '-');
    numbers = str2double(line{2});
    min = str2double(tmp1{1});
    max = str2double(tmp1{2});
    for i = min:max
        inage(i+1) = numbers/(max-min);
    end
end
fclose(f);

ingroup(1) = sum(inage(1:7)); 
ingroup(2) = sum(inage(8:19));
ingroup(3) = sum(inage(19:66));
ingroup(4) = sum(inage(66:end));


load('../AgeMixing/input_vals_group_mixing_nucdiff.mat')


weight(1) = ingroup(1)/ sum(MetaData.group==1);
weight(2) = ingroup(2)/sum(MetaData.group==2);
% sampled_adults = sum(MetaData.group==3) + sum(MetaData.group==4) + sum(MetaData.group==5);
weight(3) = ingroup(3)/sum(MetaData.group==3);
weight(4) = ingroup(3)/sum(MetaData.group==4);
weight(5) = ingroup(3)/sum(MetaData.group==5);
weight(6) = ingroup(4)/sum(MetaData.group==6);

fprintf('%.9f ', weight/sum(weight))
fprintf('\n');