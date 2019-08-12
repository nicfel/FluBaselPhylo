function [tr,boot]  = phytreeread(filename)
%PHYTREEREAD reads a NEWICK tree formatted file.
%  Nicola: now reads tree with several traits
%
%
%  TREE = PHYTREEREAD(FILENAME) reads a NEWICK tree formatted file
%  FILENAME, returning the data in the file as a PHYTREE object. FILENAME
%  can also be a URL or MATLAB character array that contains the text of a
%  NEWICK format file.
%
%  [TREE,BOOT] = PHYTREEREAD(FILENAME) returns bootstrap values when they
%  are specified either with square brackets ([]) after the branch or leaf
%  lengths or if they appear instead of the branch labels. Bootstrap values
%  that do not appear in the file default to NaN. 
%
%  The NEWICK tree format is found at:
%        http://evolution.genetics.washington.edu/phylip/newicktree.html
%
%  Note: This implementation only allows binary trees, non-binary trees
%  will be translated into a binary tree with extra branches of length 0.
%
%   Example:
%
%      tr = phytreeread('pf00002.tree')
%
%   See also GETHMMTREE, PHYTREE, PHYTREEVIEWER, PHYTREEWRITE.

% Copyright 2003-2012 The MathWorks, Inc.



if nargin==0
    [filename, pathname] = uigetfile({'*.tree';'*.dnd'},'Select Phylogenetic Tree File');
    if ~filename
        disp('Canceled, file not read.');
        tr=[];
        return;
    end
    filename = [pathname, filename];
end

% check input is char
% in a future version we may accept also cells
if ~ischar(filename)
    error(message('bioinfo:phytreeread:InvalidInput'))
end

if size(filename,1)>1  % is padded string ?
    strin = cellstr(filename);
    strin = [strin{:}];
elseif (strfind(filename(1:min(10,end)), '://')) % is an url ?
    if (~usejava('jvm'))
        error(message('bioinfo:phytreeread:NoJava'))
    end
    try
        strin = urlread(filename);
    catch allExceptions
        error(message('bioinfo:phytreeread:CannotReadURL', filename));
    end
    strin = textscan(strin,'%c','delimiter','\n');
    strin = strin{1}';
elseif  (exist(filename,'file') || ...
        exist(fullfile(pwd,filename),'file') )    %  is a valid filename ?
    fid = fopen(filename,'r');
    strin = textscan(fid,'%c','delimiter','\n');
    strin = strin{1}';
    fclose(fid);
else  % must be single a string with '\n'
    strin = textscan(filename,'%c','delimiter','\n');
    strin = strin{1}';
end

% Find quoted strings and escape control characters (G1030211)
[ins,outs]=regexp(strin,'\''[^\'']+\''','match','split');
if numel(ins)>0
    ins = regexprep(regexprep(ins,'\''',''),'([\(\):;,])','\\$1');
    strin = [outs;[ins {''}]];
    strin = sprintf('%s',strin{:});
end

% Remove inserted escape characters: 
% e.g. monkey(Africa) may be a label and it appears in the source as
% monkey\(Africa\), so these parentheses are not considered tree branches 
strinESC = regexprep(strin,'\\[\(\):;,]','__');

% temporaly mask every comma in []
mask_on=false;
for j = 1:length(strinESC)
    if strcmp(strinESC(j),'[')
        mask_on=true;
    elseif strcmp(strinESC(j),']')
        mask_on=false;
    elseif mask_on && strcmp(strinESC(j),',')
        strinESC(j)='&';
    end
end


% Characterizing the string:
numBranches = sum(strinESC==',');
numLeaves   = numBranches + 1;
numLabels   = numBranches + numLeaves;

if (numBranches == 0)
    error(message('bioinfo:phytreeread:NoCommaInInputString'))
end

% Find the string features: open and close parentheses and leaves
% e.g. strFeatures -> ((rb)((ss)((mc)w))d)
leafPositions = regexp(strinESC,'[(,][^(,)]')+1;
parenthesisPositions = regexp(strinESC,'[()]');
strFeatures = strinESC(sort([leafPositions parenthesisPositions]));


% Some consistency checking on the parenthesis
temp = cumsum((strFeatures=='(') - (strFeatures==')'));
if any(temp(1:end-1)<1) || (temp(end)~=0)
    error(message('bioinfo:phytreeread:InconsistentParentheses'))
end

dist = zeros(numLabels,1);             % allocating space for distances
tree = zeros(numBranches,2);           % allocating space for tree pointers
names = cell(numLabels,1);             % allocating space for tree labels
boot = nan(numLabels,1);             % allocating space for bootstrapping vals

% try
    
    % extract label information for the leaves
    st = regexp(strinESC,'[(,][^(,);]+','start');
    en = regexp(strinESC,'[(,][^(,);]+','end');
    for j=1:numel(st)
        elestr = strinESC(st(j)+1:en(j));
        oristr = strin(st(j)+1:en(j));
        % find out if there is a bootstrap value within []
        bootstr = regexp(elestr,'\[[\d\.]*\]','match','once');
        if ~isempty(bootstr)
            boot(j) = strread(bootstr(2:end-1),'%f');
            elestr = strrep(elestr,bootstr,'');
            oristr = strrep(oristr,bootstr,'');
        end
        coi = find(elestr==':',1,'last');        
        
        if isempty(coi) % if no colon no length, the whole label is the name
            dist(j) = 0;            
            names{j} = strin(st(j)+1:en(j));
        else % if there is colon, get name and length
            dist(j) = strread(oristr(coi+1:end),'%f');
            names{j} = oristr(1:coi-1);
        end
    end
    % uniformizing empty cells, value inside the brackets can never be empty
    % because branch names will always be empty
    [names{cellfun('isempty',names)}] = deal('');

    % extract label information for the parenthesis
    st = regexp(strinESC,')[^(,);]*','start');
    en = regexp(strinESC,')[^(,);]*','end');
    parenthesisDist = zeros(numel(st),1);
    parenthesisData = cell(numel(st),1);
    parenthesisBoot = nan(numel(st),1);
    for j=1:numel(st)
        elestr = strinESC(st(j)+1:en(j));
        oristr = strin(st(j)+1:en(j));
        % find out if there is a bootstrap value within []
        bootstr = regexp(elestr,'\[[\d\.]*\]','match','once');
        if ~isempty(bootstr)
            parenthesisBoot(j) = strread(bootstr(2:end-1),'%f');
            elestr = strrep(elestr,bootstr,'');
            oristr = strrep(oristr,bootstr,'');
        end
        coi = find(elestr==':',1,'last');
        if isempty(coi) % if no colon no length, the whole label is the name
            parenthesisDist(j) = 0;
            parenthesisData{j} = oristr;
        else % if there is colon, get name and length
            parenthesisDist(j) = strread(oristr(coi+1:end),'%f');
            parenthesisData{j} = oristr(1:coi-1);
        end
    end
    % uniformizing empty cells, value inside brackes may be empty
    if any(cellfun('isempty',parenthesisData))
        [parenthesisData{cellfun('isempty',parenthesisData)}] = deal('');
    end

    li = 1; bi = 1; pi = 1;          % indexes for leaf, branch and parentheses
    queue = zeros(1,2*numLeaves); qp = 0; % setting the queue (worst case size)

    j = 1;

    while j <= numel(strFeatures)
        switch strFeatures(j)
            case ')' % close parenthesis, pull values from the queue to create
                % a new branch and push the new branch # into the queue
                lastOpenPar = find(queue(1:qp)==0,1,'last');
                numElemInPar = min(3,qp-lastOpenPar);
                switch numElemInPar
                    case 2  % 99% of the cases, two elements in the parenthesis
                        bp = bi + numLeaves;
                        names{bp} = parenthesisData{pi};      % set name
                        dist(bp) = parenthesisDist(pi);       % set length
                        boot(bp) = parenthesisBoot(pi);       % set bootstrap
                        tree(bi,:) = queue(qp-1:qp);
                        qp = qp - 2; % writes over the open par mark
                        queue(qp) = bp;
                        bi = bi + 1;
                        pi = pi + 1;
                    case 3  % find in non-binary trees, create a phantom branch
                        bp = bi + numLeaves;
                        names{bp} = '';      % set name
                        dist(bp) = 0;        % set length
                        boot(bp) = NaN;        % set bootstrap
                        tree(bi,:) = queue(qp-1:qp);
                        qp = qp - 1; % writes over the left element
                        queue(qp) = bp;
                        bi = bi + 1;
                        j = j - 1; %repeat this closing branch to get the rest
                    case 1  % parenthesis with no meaning (holds one element)
                        qp = qp - 1;
                        queue(qp) = queue(qp+1);
                        pi = pi + 1;
                    case 0  % an empty parenthesis pair
                        error(message('bioinfo:phytreeread:ParenthesisPairWithNoData'))
                end % switch numElemInPar

            case '(' % an open parenthesis marker (0) pushed into the queue
                qp = qp + 1;
                queue(qp) = 0;

            otherwise % a new leaf pushed into the queue
                qp = qp + 1;
                queue(qp) = li;
                li = li + 1;
        end % switch strFeatures
        j = j + 1;
    end % while j ...

% catch le
%     if strcmp(le.identifier,'bioinfo:phytreeread:ParenthesisPairWithNoData')
%         rethrow(le)
%     else
%         error(message('bioinfo:phytreeread:IncorrectString'))
%     end
% end

% Remove all escape characters
names = regexprep(names,'\\([\(\):;,])','$1');

% some formats store bootstrap values instead of for branches instead of branch names,
% check here that all branch names can be converted to numbers and if so,
% consider them as bootstrap values:
if isequal(regexp(names(numLeaves+1:end),'[\d\.]+','match','once'),names(numLeaves+1:end))
    for j = numLeaves+(1:numBranches)
        if ~isempty(names{j})
           boot(j) = strread(names{j},'%f');
           names{j} = '';
        end
    end
end


% make sure all dists are greater than 0
dist = max(0,dist);

if sum(dist) == 0  % there was no distance information so force to an unitary ultrametric tree
    tr = phytree(tree,names);
elseif sum(dist(1:numLeaves)) == 0 % no dist infor for leaves, so force an ultrametric tree
    tr = phytree(tree,dist(numLeaves+1:end),names);
else % put all info into output object
    tr = phytree(tree,dist,names);
end
