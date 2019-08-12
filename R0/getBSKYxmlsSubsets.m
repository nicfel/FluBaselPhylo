% set up the bdsky xmls
clear
system('rm -r bdskysubset');
system('mkdir bdskysubset');

system('rm -r bdskysubset/xmls');
system('mkdir bdskysubset/xmls');

rep_count = 0;
for rep = 1 : 10
    for subset = 1 : 10
        clearvars -except subset rep rep_count
        % define all segments
        segments = {'HA' 'M' 'NA' 'NP' 'NS' 'PA' 'PB1' 'PB2'};

        % choose the sub-sampling probability between 20 and 90%
        sampprob = rand(1)*0.5+0.5;

        % read in the clusters
        f = fopen(['localClustersAssignment/local_cluster_sets_nucdiff_nrrep' num2str(subset) '.txt']);
        lc_nr = 1;
        nr_sams = 0;
        
        rep_count = rep_count+1

        count = 1;
        heights = zeros(0,0);
        while ~feof(f)
            line = strsplit(fgets(f),'|');
            members = strsplit(line{1},',');
            nr_sams = nr_sams + length(members);
            if length(members)>1
                for j = length(members):-1:1
                    if rand>=sampprob
                        members(j) = [];
                    end
                end
                if length(members)>1                
                    cluster_names{lc_nr,1} = ['lc_' num2str(lc_nr)];      
                    cluster_members{lc_nr,1} = members;
                    origin_height = str2double(line{end});
                    height = str2double(line{end});

                    cluster_treeHeight{lc_nr,1} = num2str(min(height,0.45));
                    cluster_origin{lc_nr,1} = num2str(min(origin_height,0.5));
                    cluster_prob{lc_nr,1} = line{2};
                    lc_nr = lc_nr+1;
                    count = count+1;
                end
            elseif length(members)==1 && rand<=sampprob
                cluster_names{lc_nr,1} = ['lc_' num2str(lc_nr)];      
                cluster_members{lc_nr,1} = members;
                origin_height = str2double(line{end});
                height = str2double(line{end});

                cluster_treeHeight{lc_nr,1} = num2str(min(height,0.45));
                cluster_origin{lc_nr,1} = num2str(min(origin_height,0.5));
                cluster_prob{lc_nr,1} = line{2};
                lc_nr = lc_nr+1;
                count = count+1;
            end
        end
        
        % count the number of sequences still there
        nr_samples = 0;
        for j = 1 : length(cluster_members)
            nr_samples = nr_samples + length(cluster_members{j});
        end
        
        % read in all segments
        for i = 1 : length(segments)
            seqs{i} = fastaread(['../Clusters/clusterNucleotidesDist/' segments{i} '.fasta']);
            bas_name{i} = cell(0,0);
            for j = 1 : length(seqs{i})
                if ~isempty(strfind(seqs{i}(j).Header,'Basel'));
                    tmp = strsplit(seqs{i}(j).Header,'|');
                    tmp = strsplit(tmp{1}, '/');
                    bas_name{i}{j} = tmp{end};
                else
                    bas_name{i}{j} = 'lala';
                end
            end
        end



        % get the sequences
        for i = 1 : length(cluster_names)
            for j = 1 : length(segments)
                for k = 1 : length(cluster_members{i})
                    ind = find(ismember(bas_name{j}, cluster_members{i}(k)));
                    Data.(cluster_names{i}).(segments{j})(k) = seqs{j}(ind);
                end
            end
        end

        % plot the cluster size distribution (sanity check)
        cl_size_for_plotting = zeros(1,100);
        for i = 1 : length(cluster_names)
            ind = length(Data.(cluster_names{i}).HA);
            cl_size_for_plotting(ind) = cl_size_for_plotting(ind)+1;
        end

        %% build the bdsky xml

        % define the cutoffs for each of the sequences
        cutoff.HA = [31 1731];
        cutoff.M = [26 1007];
        cutoff.NA = [20 1428];
        cutoff.NP = [46 1542];
        cutoff.NS = [27 864];
        cutoff.PA = [25 2175];
        cutoff.PB1 = [25 2325];
        cutoff.PB2 = [28 2307];

        % get the estimates of the mutation rates, kappa value, frequencies and 
        %TODO gamma parameters for each segment and position
        log_dat = importdata('../EvolutionaryRates/subset/combined/clusters.log');

        for j = 1 : length(segments)
            clear ind; ind = find(ismember(log_dat.textdata, ['mutationRate.' segments{j} '_1']));
            mut_1(j) = mean(log_dat.data(:,ind));
            clear ind; ind = find(ismember(log_dat.textdata, ['mutationRate.' segments{j} '_3']));
            mut_3(j) = mean(log_dat.data(:,ind));

            clear ind; ind = find(ismember(log_dat.textdata, ['kappa.' segments{j} '_1']));
            k_1(j) = mean(log_dat.data(:,ind));
            clear ind; ind = find(ismember(log_dat.textdata, ['kappa.' segments{j} '_3']));
            k_3(j) = mean(log_dat.data(:,ind));

            clear ind; ind = find(ismember(log_dat.textdata, ['gammaShape.' segments{j} '_1']));
            gam_1(j) = mean(log_dat.data(:,ind));
            clear ind; ind = find(ismember(log_dat.textdata, ['gammaShape.' segments{j} '_3']));
            gam_3(j) = mean(log_dat.data(:,ind));

            clear ind; ind = find(ismember(log_dat.textdata, ['freqParameter.' segments{j} '_11']));
            f_1(j,:) = round(mean(log_dat.data(:,ind:ind+3))/sum(mean(log_dat.data(:,ind:ind+3))),3);
            % ensure summing up to 1
            f_1(j,end) = 1 - sum(f_1(j,1:end-1));
            clear ind; ind = find(ismember(log_dat.textdata, ['freqParameter.' segments{j} '_31']));
            f_3(j,:) = round(mean(log_dat.data(:,ind:ind+3))/sum(mean(log_dat.data(:,ind:ind+3))),3);
            f_3(j,end) = 1 - sum(f_3(j,1:end-1));
        end

        clear ind; ind = find(ismember(log_dat.textdata, 'clockRate.c'));
        clock_rate = mean(log_dat.data(:,ind));    

        location = cell(0,0);
        f = fopen(['bdskysubset/bdsky_rand' num2str(rep_count) '_nr' num2str(nr_samples)  '.xml'],'w');
        fprintf(f,'<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate=''Standard'' beautistatus='''' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood:beast.mascot.dynamics:beast.mascot.distribution:beast.mascot.logger:bdsky" version="2.0">\n');

        for i = 1 : length(cluster_names)
            for j = 1 : length(segments)
                fprintf(f, '\t<data id="%s.%s">\n',cluster_names{i},segments{j});
                for k = 1 : length(Data.(cluster_names{i}).(segments{j}))
                    name = strsplit(Data.(cluster_names{i}).(segments{j})(k).Header,'|');
                    tmp2 = strsplit(strtrim(name{4}),'-');
                    if length(tmp2)==3
                        this_sequence = upper(strrep(Data.(cluster_names{i}).(segments{j})(k).Sequence(cutoff.(segments{j})(1):cutoff.(segments{j})(2)),'-','?'));
                        fprintf(f, '\t\t<sequence id="seq_%s.%s.%s" taxon="%s.%s" totalcount="4" value="%s"/>\n',...
                            cluster_names{i}, segments{j},...
                            name{1},...
                            name{1}, cluster_names{i},...
                            strrep(this_sequence, '-','?'));
                    else
                        error(['sequence ' Data.(cluster_names{i}).HA(j).Header ' shouldn''t be here']);
                    end
                end
                fprintf(f, '\t</data>\n');
            end
        end

        fprintf(f, '\n');
        fprintf(f, '<map name="Uniform" >beast.math.distributions.Uniform</map>\n');
        fprintf(f, '<map name="Exponential" >beast.math.distributions.Exponential</map>\n');
        fprintf(f, '<map name="LogNormal" >beast.math.distributions.LogNormalDistributionModel</map>\n');
        fprintf(f, '<map name="Normal" >beast.math.distributions.Normal</map>\n');
        fprintf(f, '<map name="Beta" >beast.math.distributions.Beta</map>\n');
        fprintf(f, '<map name="Gamma" >beast.math.distributions.Gamma</map>\n');
        fprintf(f, '<map name="LaplaceDistribution" >beast.math.distributions.LaplaceDistribution</map>\n');
        fprintf(f, '<map name="prior" >beast.math.distributions.Prior</map>\n');
        fprintf(f, '<map name="InverseGamma" >beast.math.distributions.InverseGamma</map>\n');
        fprintf(f, '<map name="OneOnX" >beast.math.distributions.OneOnX</map>\n');
        fprintf(f, '\n');




        fprintf(f, '<run id="mcmc" spec="MCMC" chainLength="50000000">\n');
        fprintf(f, '\t<state id="state" storeEvery="5000">\n');
        samptimes = cell(0,0);
        max_samp_times = zeros(0,0);
        min_samp_times = zeros(0,0);
        for i = 1 : length(cluster_names)
            samptimes{i} = zeros(0,0);
            for j = 1% : length(segments)
                fprintf(f, '\t\t<stateNode id="Tree.t:%s.%s" spec="beast.evolution.tree.Tree">\n', cluster_names{i}, segments{j});
                fprintf(f, '\t\t\t<trait id="dateTrait.t:%s.%s" spec="beast.evolution.tree.TraitSet" traitname="date" value="\n', cluster_names{i}, segments{j});
                for k = 1 : length(Data.(cluster_names{i}).(segments{j}))
                    tmp = strsplit(Data.(cluster_names{i}).(segments{j})(k).Header, '|');
                    tmp2 = strsplit(strtrim(tmp{4}),'-');
                    % make the time decimal
                    if length(tmp2)==3
                        deztime = (datenum(tmp{4},'yyyy-mm-dd')- datenum(tmp2{1},'yyyy'))...
                            /(datenum(num2str(str2double(tmp2{1})+1),'yyyy')-datenum(tmp2{1},'yyyy'))...
                            +str2double(tmp2{1});
                        dezstring = sprintf('%.8f', deztime);
                        samptimes{i}(end+1) = deztime;
                        name = strsplit(Data.(cluster_names{i}).(segments{j})(k).Header,'|');
                        printstring = [name{1} '.' cluster_names{i} '=' dezstring];
                        if k < length(Data.(cluster_names{i}).(segments{j}))
                            fprintf(f, '\t\t\t\t\t%s,\n', printstring);
                        else
                            fprintf(f, '\t\t\t\t\t%s',printstring);
                        end
                    else
                        error(['sequence ' Data.(cluster_names{i}).(segments{j})(k).Header ' shouldn''t be here']);
                    end
                end
                max_samp_times(i) = max(samptimes{i});
                min_samp_times(i) = min(samptimes{i});
                fprintf(f, '">\n');
                fprintf(f, '\t\t\t\t<taxa id="TaxonSet.%s.%s" spec="TaxonSet">\n', cluster_names{i}, segments{j});
                fprintf(f, '\t\t\t\t\t<alignment id="%s.%s.alignment_1" spec="FilteredAlignment" filter="1::3,2::3">\n', cluster_names{i}, segments{j});
                fprintf(f, '\t\t\t\t\t\t<data idref="%s.%s"/>\n', cluster_names{i}, segments{j});
                fprintf(f, '\t\t\t\t\t</alignment>\n');
                fprintf(f, '\t\t\t\t</taxa>\n');
                fprintf(f, '\t\t\t</trait>\n');


                fprintf(f, '\t\t\t<taxonset idref="TaxonSet.%s.%s"/>\n', cluster_names{i}, segments{j});
                fprintf(f, '\t\t</stateNode>\n');
            end
        end


        % define the times when to change rates
        rate_shifts_times = 0:2/365:max(max_samp_times)-min(min_samp_times)+0.25;

        % rate_shifts_times(end+1) = max(max_samp_times)-min(min_samp_times)+0.5;

        % get which sampling intervals have no samples
        samp_props = 0.01;
        % samp_props(rate_shifts_times>=(max(rate_shifts_times)-(max(max_samp_times)-min(min_samp_times)))) = 0.01;

        % samp_props = samp_props.*[rand(size(samp_props))*0.1+0.05];



        % fprintf(f,'\t\t<parameter id="origin" lower="0.0" name="stateNode" upper="0.2">0.0027</parameter>\n');

        for i = 1 : length(cluster_names)
        %     if strcmp(strtrim(cluster_origin{i}),'NaN')
        %         fprintf(f,'\t\t<parameter id="origin.t:%s" lower="0.0" name="stateNode" upper="0.5">0.45</parameter>\n',...
        %             cluster_names{i});
        %     else
                fprintf(f,'\t\t<parameter id="origin.t:%s" lower="0.0" name="stateNode" upper="0.1">0.01</parameter>\n',...
                    cluster_names{i});
        %     end
        end
        fprintf(f,'\t\t<parameter id="samplingProportion.t:HA" dimension="%d" lower="0.0" name="stateNode" upper="1.0">%s</parameter>\n', length(rate_shifts_times), sprintf('%f ', samp_props));
        fprintf(f,'\t\t<parameter id="becomeUninfectiousRate.t:HA" dimension="%d" lower="0.0" name="stateNode" upper="Infinity">91.2500</parameter>\n', length(rate_shifts_times));
        init_vals = normrnd(0,0.1,length(rate_shifts_times),1);
        init_vals = [1:length(rate_shifts_times)]*0.001;
        fprintf(f,'\t\t<parameter id="logReproductiveNumber.t:HA" dimension="%d" lower="-2" name="stateNode" upper="2">%s</parameter>\n', length(rate_shifts_times), sprintf('%f ', init_vals/sum(init_vals)));
        fprintf(f,'\t\t<parameter id="absoluteReproductiveNumber.t:HA" dimension="1" lower="0.0" name="stateNode" upper="15">1</parameter>\n');


        fprintf(f,'\t\t<parameter id="LogNormalSigma.R0" dimension="1" lower="0.0" name="stateNode" upper="15">0.05</parameter>\n');
        fprintf(f,'\t\t<parameter id="LogNormalSigma.Sampling" dimension="1" lower="0.0" name="stateNode" upper="15">.1</parameter>\n');

        fprintf(f, '\t\t<parameter id="clockRate.c" name="stateNode">%.12f</parameter>\n', clock_rate);

        fprintf(f, '\t\t<parameter id="oriMean" lower="0.0" name="stateNode" upper="1">0.01</parameter>\n');


        for j = 1 : length(segments)
            fprintf(f, '\t\t<parameter id="kappa.s:%s_1" lower="0.0" name="stateNode">%f</parameter>\n',segments{j},k_1(j));
        end

        for j = 1 : length(segments)
            fprintf(f, '\t\t<parameter id="kappa.s:%s_3" lower="0.0" name="stateNode">%f</parameter>\n',segments{j},k_3(j));
        end
        for j = 1 : length(segments)
            fprintf(f, '\t\t<parameter id="mutationRate.s:%s_1" name="stateNode">%f</parameter>\n',segments{j}, mut_1(j));
        end
        for j = 1 : length(segments)
            fprintf(f, '\t\t<parameter id="mutationRate.s:%s_3" name="stateNode">%f</parameter>\n',segments{j}, mut_3(j));
        end
        for j = 1 : length(segments)
            fprintf(f,'\t\t<parameter id="gammaShape.s:%s_1" name="stateNode">%f</parameter>\n',segments{j}, gam_1(j));
        end
        for j = 1 : length(segments)
            fprintf(f,'\t\t<parameter id="gammaShape.s:%s_3" name="stateNode">%f</parameter>\n',segments{j}, gam_3(j));
        end
        for j = 1 : length(segments)
            fprintf(f, '\t\t<parameter id="freqParameter.s:%s_1" dimension="4" lower="0.0" name="stateNode" upper="1.0">%s</parameter>\n',segments{j}, sprintf('%f ',f_1(j,:)));
        end
        for j = 1 : length(segments)
            fprintf(f, '\t\t<parameter id="freqParameter.s:%s_3" dimension="4" lower="0.0" name="stateNode" upper="1.0">%s</parameter>\n',segments{j}, sprintf('%f ',f_3(j,:)));
        end
        fprintf(f, '\t</state>\n\n');

        for i = 1 : length(cluster_names)
            for j = 1 %: length(segments)
                if strcmp(strtrim(cluster_origin{i}), 'NaN')
                    fprintf(f, '\t<init id="RandomTree.t:%s.%s" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:%s.%s" taxa="@%s.%s.alignment_1" rootHeight="%s">\n',cluster_names{i},segments{j},cluster_names{i},segments{j},cluster_names{i},segments{j},'0.4');
                else
                    fprintf(f, '\t<init id="RandomTree.t:%s.%s" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:%s.%s" taxa="@%s.%s.alignment_1" rootHeight="%s">\n',cluster_names{i},segments{j},cluster_names{i},segments{j},cluster_names{i},segments{j},cluster_treeHeight{i});
                end
                fprintf(f, '\t\t<populationModel id="ConstantPopulation0.t:%s.%s_1" spec="ConstantPopulation">\n',cluster_names{i},segments{j});
                fprintf(f, '\t\t\t<parameter id="randomPopSize.t:%s.%s_1" name="popSize">0.1</parameter>\n',cluster_names{i},segments{j});
                fprintf(f, '\t\t</populationModel>\n');
                fprintf(f, '\t</init>\n');
            end
        end



        fprintf(f, '\t<distribution id="posterior" spec="util.CompoundDistribution">\n');
        fprintf(f, '\t\t<distribution id="prior" spec="util.CompoundDistribution">\n');   
        fprintf(f,'\t\t\t\t<prior id="becomeUninfectiousRatePrior.t:HA" name="distribution" x="@becomeUninfectiousRate.t:HA">\n');
        fprintf(f,'\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.0" name="distr">\n');
        fprintf(f,'\t\t\t\t\t\t<parameter id="RealParameter.3" estimate="false" name="M">0.0</parameter>\n');
        fprintf(f,'\t\t\t\t\t\t<parameter id="RealParameter.4" estimate="false" name="S">1.0</parameter>\n');
        fprintf(f,'\t\t\t\t\t</LogNormal>\n');
        fprintf(f,'\t\t\t\t</prior>\n');
        % fprintf(f,'\t\t\t\t<prior id="originPrior.t:HA" name="distribution" x="@origin.t:HA">\n');
        % fprintf(f,'\t\t\t\t\t<OneOnX id="Uniform.3" name="distr"/>\n');
        % fprintf(f,'\t\t\t\t</prior>\n');
        % fprintf(f,'\t\t\t\t<prior id="reproductiveNumberPrior.t:HA" name="distribution" x="@logReproductiveNumber.t:HA">\n');
        % fprintf(f,'\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.1" name="distr">\n');
        % fprintf(f,'\t\t\t\t\t\t<parameter id="RealParameter.5" estimate="false" name="M">0.0</parameter>\n');
        % fprintf(f,'\t\t\t\t\t\t<parameter id="RealParameter.6" estimate="false" name="S">1.0</parameter>\n');
        % fprintf(f,'\t\t\t\t\t</LogNormal>\n');
        % fprintf(f,'\t\t\t\t</prior>\n');
        % fprintf(f,'\t\t\t\t<prior id="samplingProportionPrior.t:HA" name="distribution" x="@samplingProportion.t:HA">\n');
        % fprintf(f,'\t\t\t\t\t<Beta id="Beta.0" name="distr">\n');
        % fprintf(f,'\t\t\t\t\t\t<parameter id="RealParameter.1" estimate="false" name="alpha">1.0</parameter>\n');
        % fprintf(f,'\t\t\t\t\t\t<parameter id="RealParameter.2" estimate="false" name="beta">1.0</parameter>\n');
        % fprintf(f,'\t\t\t\t\t</Beta>\n');
        % fprintf(f,'\t\t\t\t</prior>\n');



        fprintf(f,'\t\t\t\t<prior id="OUsigma.t:HA" name="distribution" x="@LogNormalSigma.R0">\n');
        fprintf(f,'\t\t\t\t\t<Exponential id="LogNormalDistributionModel.1" name="distr" mean="1.0"/>\n');
        fprintf(f,'\t\t\t\t</prior>\n');

        fprintf(f,'\t\t\t\t<prior id="OUsigma.t:Sampling" name="distribution" x="@LogNormalSigma.Sampling">\n');
        fprintf(f,'\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.2" name="distr" M="0.0" S="4.0"/>\n');
        fprintf(f,'\t\t\t\t</prior>\n');
        % fprintf(f,'\t\t\t\t<prior id="OUnu.Sam.Prior" name="distribution" x="@logReproductiveNumber.t:HA">\n');
        % fprintf(f,'\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.2" name="distr" M="0.0" S="1.0"/>\n');
        % fprintf(f,'\t\t\t\t</prior>\n');





        % fprintf(f,'\t\t\t\t<prior id="OUsigma.t:HA" name="distribution" x="@LogNormalSigma.R0">\n');
        % fprintf(f,'\t\t\t\t\t<OneOnX id="Uniform.4" name="distr"/>\n');
        % fprintf(f,'\t\t\t\t</prior>\n');
        % fprintf(f,'\t\t\t\t<prior id="OUsigma.t:Sampling" name="distribution" x="@LogNormalSigma.Sampling">\n');
        % fprintf(f,'\t\t\t\t\t<OneOnX id="Uniform.5" name="distr"/>\n');
        % fprintf(f,'\t\t\t\t</prior>\n');


        for i = 1 : length(cluster_names)
            fprintf(f,'\t\t\t\t<prior id="OriginPrior.%s" x="@origin.t:%s" name="distribution">\n',cluster_names{i},cluster_names{i});
            fprintf(f,'\t\t\t\t\t<LogNormal id="oriDist.%s" name="distr" M="0.0192" S="2" meanInRealSpace="true"/>\n', cluster_names{i});
            fprintf(f,'\t\t\t\t</prior>\n');    
        end



        for i = 1 : length(cluster_names)    
            rate_shifts = rate_shifts_times + max_samp_times(i) - max(max_samp_times);
            timeOffset = 0;
            for j = 1 : length(rate_shifts)-1
                if rate_shifts(j)<=0 && rate_shifts(j+1)>0
                    rate_shifts(j) = 0.0;
                    timeOffset = j-1;
                end
            end

            rate_shifts(rate_shifts<0) = [];

            fprintf(f,'\t\t\t\t<distribution id="BirthDeathSkySerial.t:%s" spec="beast.evolution.speciation.BirthDeathSkylineModel" origin="@origin.t:%s" absoluteReproductiveNumber="@absoluteReproductiveNumber.t:HA" tree="@Tree.t:%s.%s" originIsRootEdge="true">\n', cluster_names{i}, cluster_names{i}, cluster_names{i},segments{1});
            fprintf(f,'\t\t\t\t\t\t\t<reverseTimeArrays id="reverse.%s" estimate="false" spec="parameter.BooleanParameter">true true true true</reverseTimeArrays>\n',cluster_names{i});    
            fprintf(f,'\t\t\t\t\t\t\t<parameter id="intervalTimes1.%s" estimate="false" name="birthRateChangeTimes">%s</parameter>\n',cluster_names{i}, sprintf('%.12f ', rate_shifts));    
            fprintf(f,'\t\t\t\t\t\t\t<parameter id="intervalTimes2.%s" estimate="false" name="deathRateChangeTimes">%s</parameter>\n',cluster_names{i}, sprintf('%.12f ', rate_shifts));    
            fprintf(f,'\t\t\t\t\t\t\t<parameter id="intervalTimes3.%s" estimate="false" name="samplingRateChangeTimes">%s</parameter>\n',cluster_names{i}, sprintf('%.12f ', rate_shifts));        
            fprintf(f,'\t\t\t\t\t\t\t<truncatedRealParameter name="logReproductiveNumber" spec="beast.evolution.speciation.TruncatedRealParameter" id="truncatedParam1.%s" parameter="@logReproductiveNumber.t:HA" offset="%d"/>\n',cluster_names{i}, timeOffset);    
            fprintf(f,'\t\t\t\t\t\t\t<truncatedRealParameter name="becomeUninfectiousRate" spec="beast.evolution.speciation.TruncatedRealParameter" id="truncatedParam2.%s" parameter="@becomeUninfectiousRate.t:HA" offset="%d"/>\n',cluster_names{i}, timeOffset);    
            fprintf(f,'\t\t\t\t\t\t\t<truncatedRealParameter name="samplingProportion" spec="beast.evolution.speciation.TruncatedRealParameter" id="truncatedParam3.%s" parameter="@samplingProportion.t:HA" offset="%d"/>\n',cluster_names{i}, timeOffset);    

            fprintf(f,'\t\t\t\t</distribution>\n');



        end


        [~,mrsp] = max(max_samp_times);


        fprintf(f,'\t\t\t\t<prior id="absPrior" name="distribution" x="@absoluteReproductiveNumber.t:HA">\n');
        fprintf(f,'\t\t\t\t\t<Exponential id="Uniform.4" name="distr" mean="1.0"/>\n');
        fprintf(f,'\t\t\t\t</prior>\n');


        fprintf(f,'\t\t\t\t<distribution id=''logNormal'' spec="LogNormalPrior" x="@logReproductiveNumber.t:HA" sigma="@LogNormalSigma.R0">\n');
        fprintf(f,'\t\t\t\t\t<Normal id="Uniforfsdm.4" name="x0Prior" mean="0" sigma="0.1"/>\n');
        fprintf(f,'\t\t\t\t</distribution>\n');


        fprintf(f,'\t\t\t\t<distribution spec="NonZeroPrior" x="@samplingProportion.t:HA">\n');
        fprintf(f,'\t\t\t\t\t<Beta id="Beta.0" name="x0Prior">\n');
        fprintf(f,'\t\t\t\t\t\t<parameter id="RealParameter.1" estimate="false" name="alpha">1.0</parameter>\n');
        fprintf(f,'\t\t\t\t\t\t<parameter id="RealParameter.2" estimate="false" name="beta">1.0</parameter>\n');
        fprintf(f,'\t\t\t\t\t</Beta>\n');
        fprintf(f,'\t\t\t\t</distribution>\n');


        fprintf(f,'\t\t</distribution>\n');

        fprintf(f, '\t\t<distribution id="likelihood" spec="util.CompoundDistribution" useThreads="true">\n');
        isfirst = true;
        for i = 1 : length(cluster_names)
            for j = 1 : length(segments)
                if length(cluster_members{i})>1
                    if ~isfirst
                        fprintf(f,'\t\t\t<distribution id="treeLikelihood.%s.%s._1" spec="ThreadedTreeLikelihood" branchRateModel="@StrictClock.c" tree="@Tree.t:%s.%s" useAmbiguities="true">\n',...
                            cluster_names{i},segments{j},cluster_names{i},segments{1});
                    else
                        fprintf(f,'\t\t\t<distribution id="treeLikelihood.%s.%s._1" spec="ThreadedTreeLikelihood" tree="@Tree.t:%s.%s" useAmbiguities="true">\n',...
                            cluster_names{i},segments{j},cluster_names{i},segments{1});
                    end
                    fprintf(f,'\t\t\t\t<data id="%s.%s_1" spec="FilteredAlignment" filter="1::3,2::3">\n',cluster_names{i},segments{j});
                    fprintf(f,'\t\t\t\t\t<data idref="%s.%s"/>\n',cluster_names{i},segments{j});
                    fprintf(f,'\t\t\t\t</data>\n');
                    fprintf(f,'\t\t\t\t<siteModel id="SiteModel.s:%s.%s_1" spec="SiteModel" gammaCategoryCount="4" shape="@gammaShape.s:%s_1"  mutationRate="@mutationRate.s:%s_1">\n',cluster_names{i},segments{j},segments{j},segments{j});
                    fprintf(f,'\t\t\t\t\t<parameter id="proportionInvariant.s:%s.%s_1" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>\n',cluster_names{i},segments{j});
                    fprintf(f,'\t\t\t\t\t<substModel id="hky.s:%s.%s_1" spec="HKY" kappa="@kappa.s:%s_1">\n',cluster_names{i},segments{j},segments{j});
                    fprintf(f,'\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:%s.%s_1" spec="Frequencies" frequencies="@freqParameter.s:%s_1"/>\n',cluster_names{i},segments{j},segments{j});
                    fprintf(f,'\t\t\t\t\t</substModel>\n');
                    fprintf(f,'\t\t\t\t</siteModel>\n');
                    if isfirst
                         isfirst = false;
                         fprintf(f,'\t\t\t\t<branchRateModel id="StrictClock.c" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c"/>\n');
                    end

                    fprintf(f,'\t\t\t</distribution>\n');
                    fprintf(f,'\t\t\t<distribution id="treeLikelihood.%s.%s._3" spec="ThreadedTreeLikelihood" useAmbiguities="true" branchRateModel="@StrictClock.c" tree="@Tree.t:%s.%s">\n',...
                        cluster_names{i},segments{j},cluster_names{i},segments{1});
                    fprintf(f,'\t\t\t\t<data id="%s.%s_3" spec="FilteredAlignment" filter="3::3">\n',cluster_names{i},segments{j});
                    fprintf(f,'\t\t\t\t\t<data idref="%s.%s"/>\n',cluster_names{i},segments{j});
                    fprintf(f,'\t\t\t\t</data>\n');
                    fprintf(f,'\t\t\t\t<siteModel id="SiteModel.s:%s.%s_3" spec="SiteModel"  gammaCategoryCount="4" shape="@gammaShape.s:%s_3"  mutationRate="@mutationRate.s:%s_3">\n',cluster_names{i},segments{j},segments{j},segments{j});
                    fprintf(f,'\t\t\t\t\t<parameter id="proportionInvariant.s:%s.%s_3" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>\n',cluster_names{i},segments{j});
                    fprintf(f,'\t\t\t\t\t<substModel id="hky.s:%s.%s_3" spec="HKY" kappa="@kappa.s:%s_3">\n',cluster_names{i},segments{j},segments{j});
                    fprintf(f,'\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:%s.%s_3" spec="Frequencies" frequencies="@freqParameter.s:%s_3"/>\n',cluster_names{i},segments{j},segments{j});
                    fprintf(f,'\t\t\t\t\t</substModel>\n');
                    fprintf(f,'\t\t\t\t</siteModel>\n');
                    fprintf(f,'\t\t\t</distribution>\n');
                end
            end
        end




        fprintf(f,'\t\t</distribution>\n');
        fprintf(f,'\t</distribution>\n');


        for i = 1 : length(cluster_names)
            for j = 1 %: length(segments)
                if length(cluster_members{i})>4
                    fprintf(f,'\t<operator id="CoalescentConstantTreeScaler.t:%s.%s" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:%s.%s" weight="3.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
                    fprintf(f,'\t<operator id="CoalescentConstantTreeRootScaler.t:%s.%s" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:%s.%s" weight="3.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
                    fprintf(f,'\t<operator id="CoalescentConstantUniformOperator.t:%s.%s" spec="Uniform" tree="@Tree.t:%s.%s" weight="30.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
                    fprintf(f,'\t<operator id="CoalescentConstantSubtreeSlide.t:%s.%s" spec="SubtreeSlide" tree="@Tree.t:%s.%s" weight="15.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
                    fprintf(f,'\t<operator id="CoalescentConstantNarrow.t:%s.%s" spec="Exchange" tree="@Tree.t:%s.%s" weight="15.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
                    fprintf(f,'\t<operator id="CoalescentConstantWide.t:%s.%s" spec="Exchange" isNarrow="false" tree="@Tree.t:%s.%s" weight="3.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
                    fprintf(f,'\t<operator id="CoalescentConstantWilsonBalding.t:%s.%s" spec="WilsonBalding" tree="@Tree.t:%s.%s" weight="3.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
                elseif length(cluster_members{i})>2
                    fprintf(f,'\t<operator id="CoalescentConstantTreeRootScaler.t:%s.%s" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:%s.%s" weight="3.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
                    fprintf(f,'\t<operator id="CoalescentConstantSubtreeSlide.t:%s.%s" spec="SubtreeSlide" tree="@Tree.t:%s.%s" weight="15.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
                    fprintf(f,'\t<operator id="CoalescentConstantNarrow.t:%s.%s" spec="Exchange" tree="@Tree.t:%s.%s" weight="15.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
                elseif length(cluster_members{i})>1
                    fprintf(f,'\t<operator id="CoalescentConstantTreeRootScaler.t:%s.%s" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:%s.%s" weight="3.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
                end
            end
        end




        % fprintf(f,'\t\t<operator id="OUmean.SamScaler.t:HA" spec="ScaleOperator" parameter="@OUmean.Sam" scaleFactor="0.75" weight="0.5"/>\n');
        % fprintf(f,'\t\t<operator id="OUsigma.SamScaler.t:HA" spec="ScaleOperator" parameter="@OUsigma.Sam" scaleFactor="0.75" weight="0.5"/>\n');
        % fprintf(f,'\t\t<operator id="OUnu.SamScaler.t:HA" spec="ScaleOperator" parameter="@OUnu.Sam" scaleFactor="0.75" weight="0.5"/>\n');


        fprintf(f,'\t\t<operator id="LogNormalSigma.R0Scaler.t:HA" spec="ScaleOperator" parameter="@LogNormalSigma.R0" scaleFactor="0.75" weight="5"/>\n');
        % fprintf(f,'\t\t<operator id="LogNormalSigma.SamplingScaler.t:HA" spec="ScaleOperator" parameter="@LogNormalSigma.Sampling" scaleFactor="0.75" weight="5"/>\n');

        % fprintf(f,'\t\t<operator id="OUnu.R0Scaler.t:HA" spec="ScaleOperator" parameter="@OUnu.R0" scaleFactor="0.75" weight="0.5"/>\n');
        % 
        % 
        % fprintf(f,'\t\t<operator id="samplingScaler.t:HA" spec="ScaleOperator" scaleAll="true" scaleAllIndependently="false" parameter="@samplingProportion.t:HA" scaleFactor="0.75" weight="2.0"/>\n');
        % fprintf(f,'\t\t<operator id="reproductiveNumberScaler.t:HA" spec="ScaleOperator" scaleAll="true" scaleAllIndependently="true" parameter="@logReproductiveNumber.t:HA" scaleFactor="0.75" weight="200.0"/>\n');
        % fprintf(f,'\t\t<operator id="samplingScaler.t:HA" spec="ScaleOperator" scaleAll="true" scaleAllIndependently="true" parameter="@samplingProportion.t:HA" scaleFactor="0.75" weight="200.0"/>\n');
        % fprintf(f,'\t\t<operator id="reproductiveNumberScaler.t:HA" spec="DeltaExchangeOperator" parameter="@logReproductiveNumber.t:HA" delta="0.1" weight="1000.0"/>\n');
        fprintf(f,'\t\t<operator id="samplingScaler.t:HA" spec="ScaleOperator" scaleAll="true" scaleAllIndependently="false"  parameter="@samplingProportion.t:HA" scaleFactor="0.75" weight="20.0"/>\n');

        fprintf(f,'\t\t<operator id="reproductiveNumberScaler.t:HA1" scaleNeightbours="1" spec="MultiRealRandomWalkOperator" parameter="@logReproductiveNumber.t:HA" windowSize="0.1" weight="300.0"/>\n');
        fprintf(f,'\t\t<operator id="reproductiveNumberScaler.t:HA2" scaleNeightbours="5" spec="MultiRealRandomWalkOperator" parameter="@logReproductiveNumber.t:HA" windowSize="0.1" weight="300.0"/>\n');
        fprintf(f,'\t\t<operator id="reproductiveNumberScaler.t:HA3" scaleNeightbours="20" spec="MultiRealRandomWalkOperator" parameter="@logReproductiveNumber.t:HA" windowSize="0.1" weight="300.0"/>\n');
        fprintf(f,'\t\t<operator id="reproductiveNumberScaler.t:HA4" scaleNeightbours="50" spec="MultiRealRandomWalkOperator" parameter="@logReproductiveNumber.t:HA" windowSize="0.1" weight="300.0"/>\n');


        % fprintf(f,'\t\t<operator id="absReproductiveNumberScaler.t:HA" spec="ScaleOperator" parameter="@absoluteReproductiveNumber.t:HA" scaleFactor="0.75" weight="20.0"/>\n');

        % fprintf(f,'\t\t<operator id="oriMeanScaler" spec="ScaleOperator" parameter="@oriMean" scaleFactor="0.75" weight="5.0"/>\n');

        % fprintf(f,'\t\t<operator id="origScaler" spec="ScaleOperator" parameter="@origin" scaleFactor="0.75" weight="20.0"/>\n');

        for i = 1 : length(cluster_names)
            fprintf(f,'\t\t<operator id="origScaler.t:%s" spec="ScaleOperator" parameter="@origin.t:%s" scaleFactor="0.75" weight="2.0"/>\n',cluster_names{i},cluster_names{i});
        end


%          for i = 1 : length(cluster_names)
%              if length(cluster_members{i})>1
%                 fprintf(f,'\t<logger id="treelog.%s" fileName="%s_subs%d.trees" logEvery="50000" mode="tree">\n', cluster_names{i}, cluster_names{i}, subset);
%                 fprintf(f,'\t\t<log id="logTrees.%s" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:%s.%s"/>\n',cluster_names{i},cluster_names{i},segments{1});
%                 fprintf(f,'\t</logger>\n');
%              end
%          end


        fprintf(f,'\t<logger id="tracelog.clusters" fileName="clusters_rand%d_nr%d.log" logEvery="50000" model="@posterior" sanitiseHeaders="true" sort="smart">\n', rep_count, nr_samples);
        fprintf(f,'\t\t<log idref="posterior"/>\n');
        fprintf(f,'\t\t<log idref="likelihood"/>\n');
        fprintf(f,'\t\t<log idref="prior"/>\n');
        % fprintf(f,'\t\t<log idref="origin"/>\n');


        for i = 1 : length(cluster_names)
            for j = 1 : length(segments)
                if length(cluster_members{i})>1
                    fprintf(f,'\t\t<log idref="treeLikelihood.%s.%s._1"/>\n', cluster_names{i},segments{j});
                    fprintf(f,'\t\t<log idref="treeLikelihood.%s.%s._3"/>\n', cluster_names{i},segments{j});
                end
            end
        end
        for i = 1 : length(cluster_names)
            fprintf(f,'\t\t<log idref="BirthDeathSkySerial.t:%s"/>\n', cluster_names{i});
        end
        fprintf(f,'\t</logger>\n');

        fprintf(f,'\t<logger id="tracelog.origin" fileName="heights_rand%d_nr%d.log" logEvery="50000" model="@posterior" sanitiseHeaders="true" sort="smart">\n', rep_count, nr_samples);

        fprintf(f,'\t\t<log idref="oriMean"/>\n');
        for i = 1 : length(cluster_names)
            fprintf(f,'\t\t<log idref="origin.t:%s"/>\n', cluster_names{i});
        end
        for i = 1 : length(cluster_names)
            for j = 1% : length(segments)
                if length(cluster_members{i})>1
                    fprintf(f,'\t\t<log id="TreeHeight.t:%s.%s" spec="beast.evolution.tree.TreeHeightLogger" tree="@Tree.t:%s.%s"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
                end
            end    
        end    

        fprintf(f,'\t</logger>\n');


        fprintf(f,'\t<logger id="tracelog.bdsky" fileName="bdsky_rand%d_nr%d_sub%d.log" logEvery="50000" model="@posterior" sanitiseHeaders="true" sort="smart">\n', rep_count, nr_samples, subset);
        fprintf(f,'\t\t<log idref="LogNormalSigma.R0"/>\n');
        fprintf(f,'\t\t<log idref="LogNormalSigma.Sampling"/>\n');
        fprintf(f,'\t\t<log idref="absoluteReproductiveNumber.t:HA"/>\n');
        fprintf(f,'\t\t<log idref="samplingProportion.t:HA"/>\n');
        fprintf(f,'\t\t<log idref="becomeUninfectiousRate.t:HA"/>\n');
        fprintf(f,'\t\t<log idref="logReproductiveNumber.t:HA"/>\n');
        fprintf(f,'\t</logger>\n');

        fprintf(f,'\t<logger id="screenlog" logEvery="100000">\n');
        fprintf(f,'\t\t<log idref="posterior"/>\n');
        fprintf(f,'\t\t<log id="ESS.0" spec="util.ESS" arg="@posterior"/>\n');
        fprintf(f,'\t\t<log idref="likelihood"/>\n');
        fprintf(f,'\t\t<log idref="prior"/>\n');
        fprintf(f,'\t</logger>\n');
        fprintf(f,'</run>\n');
        fprintf(f,'</beast>\n');

        fclose(f);
        disp('done');
    end
end
    
%% make replicates
fclose('all')

xmls = dir('bdskysubset/*.xml');

for i = 1 : length(xmls)
    for r = 0 : 2

        f = fopen(['bdskysubset/' xmls(i).name]);
        g = fopen(['bdskysubset/xmls/' strrep(xmls(i).name, '.xml', sprintf('_rep%d.xml', r))], 'w');

        while ~feof(f)
            line = fgets(f);
            if ~isempty(strfind(line, 'fileName="'))
                fprintf(g, strrep(line, 'fileName="', ['fileName="rep' num2str(r) '_']));
            else
                fprintf(g, '%s', line);
            end
        end
        fclose('all')
    end
end

