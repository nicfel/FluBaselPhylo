function [] = getNucleotideDifference(use_basel_ind)


if exist('input_distance.mat')>0
    load('input_distance.mat')
else
    % read in the sequences
    fastas = dir('combined/*.afasta');
    for i = 1 : length(fastas)
        if length(fastas(i).name) < 12
            Data{i} = fastaread(['combined/' fastas(i).name]);
        end
    end
    % make a "library" with all the headers
    for i = 1 : length(Data)
        for j = 1 : length(Data{i})
            id{i}{j} = Data{i}(j).Header;
        end
    end
    % get all Basel id's
    c = 1;
    for j = 1 : length(id{1})
        if ~isempty(strfind(id{1}{j},'Basel'))
            basel_ids{c} = id{1}{j};
            c = c+1;
        end
    end
    uni_id = cell(1,0);
    for i=1:length(id)
        uni_id = unique([uni_id id{i}]);
    end
    save('input_distance.mat')
end
%%
if exist('distances.mat')>0
    load('distances.mat')
else
    baselDist_tmp = zeros(1,0);
    a = use_basel_ind;
    disp(a)
    ind_a = zeros(1,length(Data));
    for i = 1 : length(Data)
        ind = find(ismember(id{i}, basel_ids{a}));
        if ~isempty(ind)
            ind_a(i) = ind;            
        end
    end

    for b = 1 : length(uni_id)
        ind_b = zeros(1,length(Data));
        for i = 1 : length(Data)
            ind = find(ismember(id{i}, uni_id{b}));
            if ~isempty(ind)
                ind_b(i) = ind;
            end
        end

        diff = 0;
        all_length = 0;
        for i = 1 : length(Data)
            if ind_a(i)>0 & ind_b(i)>0
                seq1 = Data{i}(ind_a(i)).Sequence;
                seq2 = Data{i}(ind_b(i)).Sequence;
                seq1 = strrep(seq1, '-', 'N');
                seq2 = strrep(seq2, '-', 'N');


                differences = seq1~=seq2;
                allIsN = seq1=='N' | seq2=='N';
                isN = seq1(differences)=='N' | seq2(differences)=='N';
                diff = diff + sum(isN==0);
                all_length = all_length + sum(allIsN==0);
            end
        end
        baselDist_tmp(b) = diff/all_length;
        baselDist_tmp(b) = diff/all_length;
    end
    f = fopen(['dist_' num2str(use_basel_ind) '.csv'], 'w');
    for i = 1 : length(baselDist_tmp)
        fprintf(f, '%f\n', baselDist_tmp(i));
    end
    fclose(f);
    return;

end


