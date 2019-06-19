clear 
% define the segments

segments = {'HA' 'M' 'NA' 'NP' 'NS' 'PA' 'PB1' 'PB2'};

% define the cutoffs for each of the sequences
cutoff.HA = [31 1731];
cutoff.M = [26 1007];
cutoff.NA = [20 1428];
cutoff.NP = [46 1542];
cutoff.NS = [27 864];
cutoff.PA = [25 2175];
cutoff.PB1 = [25 2325];
cutoff.PB2 = [28 2307];


%% read in the sequences again
for j = 1 : length(segments) 
    Data.(segments{j}) = fastaread(['nonaligned/' segments{j} '.fasta']);
end


%%


mut_1 = [0.55, 0.38, 0.705,	0.25, 0.92, 0.318, 0.319, 0.366];
mut_3 = [2.035, 1.244, 2.352, 2.126, 1.131, 2.244, 2.498, 2.357];

k_1 = [7.474, 8.409, 7.59, 7.515, 10.125, 7.433, 5.824, 13.777];
k_3 = [9.118,13.863,11.58,12.525,9.841,13.883,17.098,12.6];


f_1 = [0.359,0.173,0.241,0.227;
    0.295,0.225,0.262,0.218;
    0.329,0.175,0.261,0.235;
    0.331,0.187,0.279,0.203;
    0.346,0.193,0.242,0.219;
    0.35,0.177,0.233,0.24;
    0.369,0.188,0.224,0.219;
    0.335,0.182,0.254,0.229];
f_3 = [0.317,0.23,0.181,0.272;
    0.27,0.171,0.262,0.297;
    0.28,0.229,0.183,0.308;
    0.332,0.198,0.229,0.241;
    0.329,0.184,0.199,0.288;
    0.326,0.211,0.23,0.233;
    0.337,0.211,0.231,0.221;
    0.364,0.18,0.224,0.232];
    
system('rm -r subset');
system('mkdir subset');

cluster_names = {'c1'};
for i = 1 : length(cluster_names)


    f = fopen(['subset/subset.xml'],'w');
    fprintf(f,'<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate=''Standard'' beautistatus='''' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" required="" version="2.4">\n');

    for j = 1 : length(segments)
        fprintf(f, '\t<data id="%s.%s">\n',cluster_names{i},segments{j});
        for k = 1 : length(Data.(segments{j}))
            name = strsplit(Data.(segments{j})(k).Header,'|');
            tmp2 = strsplit(strtrim(name{4}),'-');
            if length(tmp2)==3
                this_sequence = upper(strrep(Data.(segments{j})(k).Sequence(cutoff.(segments{j})(1):cutoff.(segments{j})(2)),'-','?'));
                fprintf(f, '\t\t<sequence id="seq_%s.%s.%s" taxon="%s.%s" totalcount="4" value="%s"/>\n',...
                    cluster_names{i}, segments{j},...
                    name{1},...
                    name{1}, cluster_names{i},...
                    strrep(this_sequence, '-','?'));
            else
                error(['sequence ' Data.HA(j).Header ' shouldn''t be here']);
            end
        end
        fprintf(f, '\t</data>\n');
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


    fprintf(f, '<run id="mcmc" spec="MCMC" chainLength="100000000">\n');
    fprintf(f, '\t<state id="state" storeEvery="5000">\n');
    for j = 1 : length(segments)
        fprintf(f, '\t\t<tree id="Tree.t:%s.%s" name="stateNode">\n', cluster_names{i}, segments{j});
        fprintf(f, '\t\t\t<trait id="dateTrait.t:%s.%s" spec="beast.evolution.tree.TraitSet" traitname="date" value="\n', cluster_names{i}, segments{j});
        for k = 1 : length(Data.(segments{j}))
            tmp = strsplit(Data.(segments{j})(k).Header, '|');
            tmp2 = strsplit(strtrim(tmp{4}),'-');
            % make the time decimal
            if length(tmp2)==3
                deztime = (datenum(tmp{4},'yyyy-mm-dd')- datenum(tmp2{1},'yyyy'))...
                    /(datenum(num2str(str2double(tmp2{1})+1),'yyyy')-datenum(tmp2{1},'yyyy'))...
                    +str2double(tmp2{1});
                dezstring = sprintf('%.8f', deztime);
                name = strsplit(Data.(segments{j})(k).Header,'|');
                printstring = [name{1} '.' cluster_names{i} '=' dezstring];
                if k < length(Data.(segments{j}))
                    fprintf(f, '\t\t\t\t\t%s,\n', printstring);
                else
                    fprintf(f, '\t\t\t\t\t%s',printstring);
                end
            else
                error(['sequence ' Data.(segments{j})(k).Header ' shouldn''t be here']);
            end
        end
        fprintf(f, '">\n');
        fprintf(f, '\t\t\t\t<taxa id="TaxonSet.%s.%s" spec="TaxonSet">\n', cluster_names{i}, segments{j});
        fprintf(f, '\t\t\t\t\t<alignment id="%s.%s.alignment_1" spec="FilteredAlignment" filter="1::3,2::3">\n', cluster_names{i}, segments{j});
        fprintf(f, '\t\t\t\t\t\t<data idref="%s.%s"/>\n', cluster_names{i}, segments{j});
        fprintf(f, '\t\t\t\t\t</alignment>\n');
        fprintf(f, '\t\t\t\t</taxa>\n');
        fprintf(f, '\t\t\t</trait>\n');
        fprintf(f, '\t\t\t<taxonset idref="TaxonSet.%s.%s"/>\n', cluster_names{i}, segments{j});
        fprintf(f, '\t\t</tree>\n');
    end


    fprintf(f, '\t\t<parameter id="clockRate.c" name="stateNode">0.0036</parameter>\n');
    fprintf(f, '\t\t<parameter id="popSize.t" name="stateNode">0.3</parameter>\n');

    for j = 1 : length(segments)
        fprintf(f, '\t\t<parameter id="kappa.s:%s_1" lower="0.0" name="stateNode">%.3f</parameter>\n',segments{j},k_1(j));
    end

    for j = 1 : length(segments)
        fprintf(f, '\t\t<parameter id="kappa.s:%s_3" lower="0.0" name="stateNode">%.3f</parameter>\n',segments{j},k_3(j));
    end
    for j = 1 : length(segments)
        fprintf(f, '\t\t<parameter id="mutationRate.s:%s_1" name="stateNode">%.3f</parameter>\n',segments{j}, mut_1(j));
    end
    for j = 1 : length(segments)
        fprintf(f, '\t\t<parameter id="mutationRate.s:%s_3" name="stateNode">%.3f</parameter>\n',segments{j}, mut_3(j));
    end
  
    for j = 1 : length(segments)
        fprintf(f,'\t\t<parameter id="gammaShape.s:%s_1" name="stateNode">1.0</parameter>\n',segments{j});
    end

    
    for j = 1 : length(segments)
        fprintf(f,'\t\t<parameter id="gammaShape.s:%s_3" name="stateNode">1.0</parameter>\n',segments{j});
    end
    
    for j = 1 : length(segments)
        fprintf(f, '\t\t<parameter id="freqParameter.s:%s_1" dimension="4" lower="0.0" name="stateNode" upper="1.0">%s</parameter>\n',segments{j}, sprintf('%.3f ',f_1(j,:)));
    end
    
    for j = 1 : length(segments)
        fprintf(f, '\t\t<parameter id="freqParameter.s:%s_3" dimension="4" lower="0.0" name="stateNode" upper="1.0">%s</parameter>\n',segments{j}, sprintf('%.3f ',f_3(j,:)));
    end
    fprintf(f, '\t</state>\n\n');

    for j = 1 : length(segments)
        fprintf(f, '\t<init id="RandomTree.t:%s.%s" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:%s.%s" taxa="@%s.%s.alignment_1">\n',cluster_names{i},segments{j},cluster_names{i},segments{j},cluster_names{i},segments{j});
        fprintf(f, '\t\t<populationModel id="ConstantPopulation0.t:%s.%s_1" spec="ConstantPopulation">\n',cluster_names{i},segments{j});
        fprintf(f, '\t\t\t<parameter id="randomPopSize.t:%s.%s_1" name="popSize">1.0</parameter>\n',cluster_names{i},segments{j});
        fprintf(f, '\t\t</populationModel>\n');
        fprintf(f, '\t</init>\n');
    end

    fprintf(f, '\t<distribution id="posterior" spec="util.CompoundDistribution">\n');
    fprintf(f, '\t\t<distribution id="prior" spec="util.CompoundDistribution">\n');
    for j = 1 : length(segments)
        fprintf(f, '\t\t\t<distribution id="CoalescentConstant.t:%s.%s" spec="Coalescent">\n',cluster_names{i},segments{j});
        fprintf(f, '\t\t\t\t<populationModel id="ConstantPopulation.t:%s.%s" spec="ConstantPopulation" popSize="@popSize.t"/>\n',cluster_names{i},segments{j});
        fprintf(f, '\t\t\t\t<treeIntervals id="TreeIntervals.t:%s.%s" spec="TreeIntervals" tree="@Tree.t:%s.%s"/>\n',cluster_names{i},segments{j},cluster_names{i},segments{j});
        fprintf(f, '\t\t\t</distribution>\n');
    end

    fprintf(f, '\t\t\t<prior id="ClockPrior.c" name="distribution" x="@clockRate.c">\n');
    fprintf(f, '\t\t\t\t<Uniform id="Uniform.0" name="distr" upper="Infinity"/>\n');
    fprintf(f, '\t\t\t</prior>\n');
    for j = 1 : length(segments)
        fprintf(f, '\t\t\t<prior id="KappaPrior.s:%s_1" name="distribution" x="@kappa.s:%s_1">\n', segments{j}, segments{j});
        fprintf(f, '\t\t\t\t<LogNormal id="LogNormalDistributionModel.%d_1" name="distr">\n', j);
        fprintf(f, '\t\t\t\t\t<parameter id="RealParameter.1.%s" estimate="false" name="M">1.0</parameter>\n',segments{j});
        fprintf(f, '\t\t\t\t\t<parameter id="RealParameter.2.%s" estimate="false" name="S">1.25</parameter>\n',segments{j});
        fprintf(f, '\t\t\t\t</LogNormal>\n');
        fprintf(f, '\t\t\t</prior>\n');
    end

    for j = 1 : length(segments)
        fprintf(f, '\t\t\t<prior id="KappaPrior.s:%s_3" name="distribution" x="@kappa.s:%s_3">\n', segments{j}, segments{j});
        fprintf(f, '\t\t\t\t<LogNormal id="LogNormalDistributionModel.%d_3" name="distr">\n', j);
        fprintf(f, '\t\t\t\t\t<parameter id="RealParameter.5.%s" estimate="false" name="M">1.0</parameter>\n',segments{j});
        fprintf(f, '\t\t\t\t\t<parameter id="RealParameter.6.%s" estimate="false" name="S">1.25</parameter>\n',segments{j});
        fprintf(f, '\t\t\t\t</LogNormal>\n');
        fprintf(f, '\t\t\t</prior>\n');
    end
    
    
    for j = 1 : length(segments)
        fprintf(f, '\t\t\t<prior id="GammaPrior.s:%s_1" name="distribution" x="@gammaShape.s:%s_1">\n', segments{j}, segments{j});
        fprintf(f, '\t\t\t\t<Exponential id="ExponentialDistribution.%d_3" name="distr">\n', j);
        fprintf(f, '\t\t\t\t\t<parameter id="RealParameter.50.%s" estimate="false" name="mean">1.0</parameter>\n',segments{j});
        fprintf(f, '\t\t\t\t</Exponential>\n');
        fprintf(f, '\t\t\t</prior>\n');
    end
    for j = 1 : length(segments)
        fprintf(f, '\t\t\t<prior id="GammaPrior.s:%s_3" name="distribution" x="@gammaShape.s:%s_3">\n', segments{j}, segments{j});
        fprintf(f, '\t\t\t\t<Exponential id="ExponentialDistribution.%d_30" name="distr">\n', j);
        fprintf(f, '\t\t\t\t\t<parameter id="RealParameter.20.%s" estimate="false" name="mean">1.0</parameter>\n',segments{j});
        fprintf(f, '\t\t\t\t</Exponential>\n');
        fprintf(f, '\t\t\t</prior>\n');
    end


    
    fprintf(f, '\t\t\t<prior id="PopSizePrior.t" name="distribution" x="@popSize.t">\n');
    fprintf(f, '\t\t\t\t<OneOnX id="OneOnX.1" name="distr"/>\n');
    fprintf(f, '\t\t\t</prior>\n');
    fprintf(f, '\t\t</distribution>\n');
    fprintf(f, '\t\t<distribution id="likelihood" spec="util.CompoundDistribution" useThreads="true">\n');
    isfirst = true;
    for j = 1 : length(segments)
        if ~isfirst
            fprintf(f,'\t\t\t<distribution id="treeLikelihood.%s.%s._1" spec="ThreadedTreeLikelihood" branchRateModel="@StrictClock.c" tree="@Tree.t:%s.%s" useAmbiguities="true">\n',...
                cluster_names{i},segments{j},cluster_names{i},segments{j});
        else
            fprintf(f,'\t\t\t<distribution id="treeLikelihood.%s.%s._1" spec="ThreadedTreeLikelihood" tree="@Tree.t:%s.%s" useAmbiguities="true">\n',...
                cluster_names{i},segments{j},cluster_names{i},segments{j});
        end
        fprintf(f,'\t\t\t\t<data id="%s.%s_1" spec="FilteredAlignment" filter="1::3,2::3">\n',cluster_names{i},segments{j});
        fprintf(f,'\t\t\t\t\t<data idref="%s.%s"/>\n',cluster_names{i},segments{j});
        fprintf(f,'\t\t\t\t</data>\n');
        fprintf(f,'\t\t\t\t<siteModel id="SiteModel.s:%s.%s_1" spec="SiteModel"  gammaCategoryCount="4" shape="@gammaShape.s:%s_1" mutationRate="@mutationRate.s:%s_1">\n',cluster_names{i},segments{j},segments{j},segments{j});
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
            cluster_names{i},segments{j},cluster_names{i},segments{j});
        fprintf(f,'\t\t\t\t<data id="%s.%s_3" spec="FilteredAlignment" filter="3::3">\n',cluster_names{i},segments{j});
        fprintf(f,'\t\t\t\t\t<data idref="%s.%s"/>\n',cluster_names{i},segments{j});
        fprintf(f,'\t\t\t\t</data>\n');
        fprintf(f,'\t\t\t\t<siteModel id="SiteModel.s:%s.%s_3" spec="SiteModel" gammaCategoryCount="4" shape="@gammaShape.s:%s_1" mutationRate="@mutationRate.s:%s_3">\n',cluster_names{i},segments{j},segments{j},segments{j});
        fprintf(f,'\t\t\t\t\t<parameter id="proportionInvariant.s:%s.%s_3" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>\n',cluster_names{i},segments{j});
        fprintf(f,'\t\t\t\t\t<substModel id="hky.s:%s.%s_3" spec="HKY" kappa="@kappa.s:%s_3">\n',cluster_names{i},segments{j},segments{j});
        fprintf(f,'\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:%s.%s_3" spec="Frequencies" frequencies="@freqParameter.s:%s_3"/>\n',cluster_names{i},segments{j},segments{j});
        fprintf(f,'\t\t\t\t\t</substModel>\n');
        fprintf(f,'\t\t\t\t</siteModel>\n');
        fprintf(f,'\t\t\t</distribution>\n');
    end

    fprintf(f,'\t\t</distribution>\n');
    fprintf(f,'\t</distribution>\n');

    fprintf(f,'\t<operator id="StrictClockRateScaler.c:HA_1,2" spec="ScaleOperator" parameter="@clockRate.c" scaleFactor="0.75" weight="3.0"/>\n');


    fprintf(f,'\t<operator id="strictClockUpDownOperator.c:HA_1,2" spec="UpDownOperator" scaleFactor="0.75" weight="3.0">\n');
    fprintf(f,'\t\t<up idref="clockRate.c"/>\n');
    for j = 1 : length(segments)
        fprintf(f,'\t\t<down idref="Tree.t:%s.%s"/>\n', cluster_names{i},segments{j});
    end

    fprintf(f,'\t</operator>\n');
    for j = 1 : length(segments)
        fprintf(f,'\t<operator id="KappaScaler.s:%s_1" spec="ScaleOperator" parameter="@kappa.s:%s_1" scaleFactor="0.5" weight="0.1"/>\n', segments{j}, segments{j});
    end

    for j = 1 : length(segments)
        fprintf(f,'\t<operator id="KappaScaler.s:%s_3" spec="ScaleOperator" parameter="@kappa.s:%s_3" scaleFactor="0.5" weight="0.1"/>\n', segments{j}, segments{j});
    end

    for j = 1 : length(segments)
        fprintf(f,'\t<operator id="FrequenciesExchanger.s:%s_1" spec="DeltaExchangeOperator" delta="0.01" weight="0.1">\n', segments{j});
        fprintf(f,'\t\t<parameter idref="freqParameter.s:%s_1"/>\n', segments{j});
        fprintf(f,'\t</operator>\n');
    end
    for j = 1 : length(segments)
        fprintf(f,'\t<operator id="FrequenciesExchanger.s:%s_3" spec="DeltaExchangeOperator" delta="0.01" weight="0.1">\n', segments{j});
        fprintf(f,'\t\t<parameter idref="freqParameter.s:%s_3"/>\n', segments{j});
        fprintf(f,'\t</operator>\n');
    end
    
    for j = 1 : length(segments)
        fprintf(f,'\t<operator id="alpha_scaler_1.%s" spec="ScaleOperator" parameter="@gammaShape.s:%s_1" scaleFactor="0.75" weight="3.0"/>\n', segments{j}, segments{j});
    end
    for j = 1 : length(segments)
        fprintf(f,'\t<operator id="alpha_scaler_3.%s" spec="ScaleOperator" parameter="@gammaShape.s:%s_3" scaleFactor="0.75" weight="3.0"/>\n', segments{j}, segments{j});
    end




    fprintf(f,'\t<operator id="FixMeanMutationRatesOperator" spec="DeltaExchangeOperator" delta="0.75" weight="2.0">\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:HA_1"/>\n');
    % fprintf(f,'\t\t<parameter idref="mutationRate.s:HA_2"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:HA_3"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:M_1"/>\n');
    % fprintf(f,'\t\t<parameter idref="mutationRate.s:M_2"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:M_3"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:NA_1"/>\n');
    % fprintf(f,'\t\t<parameter idref="mutationRate.s:NA_2"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:NA_3"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:NP_1"/>\n');
    % fprintf(f,'\t\t<parameter idref="mutationRate.s:NP_2"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:NP_3"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:NS_1"/>\n');
    % fprintf(f,'\t\t<parameter idref="mutationRate.s:NS_2"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:NS_3"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:PA_1"/>\n');
    % fprintf(f,'\t\t<parameter idref="mutationRate.s:PA_2"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:PA_3"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:PB1_1"/>\n');
    % fprintf(f,'\t\t<parameter idref="mutationRate.s:PB1_2"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:PB1_3"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:PB2_1"/>\n');
    % fprintf(f,'\t\t<parameter idref="mutationRate.s:PB2_2"/>\n');
    fprintf(f,'\t\t<parameter idref="mutationRate.s:PB2_3"/>\n');
    fprintf(f,'\t\t<weightvector id="weightparameter" spec="parameter.IntegerParameter" dimension="16" estimate="false" lower="0" upper="0">1134 567 654 327 940 470 998 499 558 279 1434 717 1516 758 1520 760</weightvector>\n');
    fprintf(f,'\t</operator>\n');

    for j = 1 : length(segments)
        fprintf(f,'\t<operator id="CoalescentConstantTreeScaler.t:%s.%s" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:%s.%s" weight="3.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
        fprintf(f,'\t<operator id="CoalescentConstantTreeRootScaler.t:%s.%s" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:%s.%s" weight="3.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
        fprintf(f,'\t<operator id="CoalescentConstantUniformOperator.t:%s.%s" spec="Uniform" tree="@Tree.t:%s.%s" weight="30.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
        fprintf(f,'\t<operator id="CoalescentConstantSubtreeSlide.t:%s.%s" spec="SubtreeSlide" tree="@Tree.t:%s.%s" weight="15.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
        fprintf(f,'\t<operator id="CoalescentConstantNarrow.t:%s.%s" spec="Exchange" tree="@Tree.t:%s.%s" weight="15.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
        fprintf(f,'\t<operator id="CoalescentConstantWide.t:%s.%s" spec="Exchange" isNarrow="false" tree="@Tree.t:%s.%s" weight="3.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
        fprintf(f,'\t<operator id="CoalescentConstantWilsonBalding.t:%s.%s" spec="WilsonBalding" tree="@Tree.t:%s.%s" weight="3.0"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
    end


    fprintf(f,'\t<operator id="PopSizeScaler.t" spec="ScaleOperator" parameter="@popSize.t" scaleFactor="0.75" weight="3.0"/>\n');


    for j = 1 : length(segments)
        fprintf(f,'\t<logger id="treelog.t:%s.%s" fileName="%s.%s.trees" logEvery="1000000" mode="tree">\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
        fprintf(f,'\t\t<log id="TreeWithMetaDataLogger.t:%s.%s" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:%s.%s"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
        fprintf(f,'\t</logger>\n');
    end



    fprintf(f,'\t<logger id="tracelog" fileName="clusters.log" logEvery="1000000" model="@posterior" sanitiseHeaders="true" sort="smart">\n');
    fprintf(f,'\t\t<log idref="posterior"/>\n');
    fprintf(f,'\t\t<log idref="likelihood"/>\n');
    fprintf(f,'\t\t<log idref="prior"/>\n');
    for j = 1 : length(segments)
        fprintf(f,'\t\t<log idref="treeLikelihood.%s.%s._1"/>\n', cluster_names{i},segments{j});
        fprintf(f,'\t\t<log idref="treeLikelihood.%s.%s._3"/>\n', cluster_names{i},segments{j});
    end


    for j = 1 : length(segments)
        fprintf(f,'\t\t<log idref="CoalescentConstant.t:%s.%s"/>\n', cluster_names{i}, segments{j});
    end


    for j = 1 : length(segments)
        fprintf(f,'\t\t<log id="TreeHeight.t:%s.%s" spec="beast.evolution.tree.TreeHeightLogger" tree="@Tree.t:%s.%s"/>\n', cluster_names{i}, segments{j}, cluster_names{i}, segments{j});
    end
    
    fprintf(f,'\t\t<log idref="clockRate.c"/>\n');
    fprintf(f,'\t\t<log idref="popSize.t"/>\n');

    for j = 1 : length(segments)
        fprintf(f,'\t\t<log idref="kappa.s:%s_1"/>\n', segments{j});
    end
    % for j = 1 : length(segments)
    %     fprintf(f,'\t\t<log idref="kappa.s:%s_2"/>\n', segments{j});
    % end
    for j = 1 : length(segments)
        fprintf(f,'\t\t<log idref="kappa.s:%s_3"/>\n', segments{j});
    end
    for j = 1 : length(segments)
        fprintf(f,'\t\t<log idref="mutationRate.s:%s_1"/>\n', segments{j});
    end
    % for j = 1 : length(segments)
    %     fprintf(f,'\t\t<log idref="mutationRate.s:%s_2"/>\n', segments{j});
    % end
    for j = 1 : length(segments)
        fprintf(f,'\t\t<log idref="mutationRate.s:%s_3"/>\n', segments{j});
    end
    
    for j = 1 : length(segments)
        fprintf(f,'\t\t<log idref="gammaShape.s:%s_1" />\n', segments{j});
    end
    for j = 1 : length(segments)
        fprintf(f,'\t\t<log idref="gammaShape.s:%s_3" />\n', segments{j});
    end


    
    for j = 1 : length(segments)
        fprintf(f,'\t\t<log idref="freqParameter.s:%s_1"/>\n', segments{j});
    end
    % for j = 1 : length(segments)
    %     fprintf(f,'\t\t<log idref="freqParameter.s:%s_2"/>\n', segments{j});
    % end
    for j = 1 : length(segments)
        fprintf(f,'\t\t<log idref="freqParameter.s:%s_3"/>\n', segments{j});
    end
    fprintf(f,'\t</logger>\n');

    fprintf(f,'\t<logger id="screenlog" logEvery="1000000">\n');
    fprintf(f,'\t\t<log idref="posterior"/>\n');
    fprintf(f,'\t\t<log id="ESS.0" spec="util.ESS" arg="@posterior"/>\n');
    fprintf(f,'\t\t<log idref="likelihood"/>\n');
    fprintf(f,'\t\t<log idref="prior"/>\n');
    fprintf(f,'\t</logger>\n');
    fprintf(f,'</run>\n');
    fprintf(f,'</beast>\n');

    fclose(f);
end

%% make replictes of the beast xmls
system('rm -r subset/xmls');
system('mkdir subset/xmls');
for i = 0 : 9
    f = fopen('subset/subset.xml');
    g = fopen(sprintf('subset/xmls/subset_rep%d.xml',i), 'w');
    while ~feof(f)
        line = fgets(f);
        if ~isempty(strfind(line, 'name="popSize">1.0</parameter>'))
            fprintf(g, strrep(line, 'name="popSize">1.0</parameter>', sprintf('name="popSize">%.5f</parameter>', exprnd(1)) ));
        elseif ~isempty(strfind(line, 'fileName="c1'))
            fprintf(g, strrep(line, 'fileName="c1.', ['fileName="rep' num2str(i) '_']));
        elseif ~isempty(strfind(line, 'fileName="clus'))
            fprintf(g, strrep(line, 'fileName="', ['fileName="rep' num2str(i) '_']));
        else
            fprintf(g, line);
        end
    end
    fclose(f);
    fclose(g);
end
system('rm -r subset/subset.xml');
