function  trials = selecttrial(TE,selSTR,varargin)
%SELECTTRIAL   Execute logical expression to filter trials.
%   TRIALS = SELECTTRIAL(TE,SELSTR) returns the trial indices (for trial
%   structure TE) satisfying the required criteria set up by SELSTR, which
%   should be a logical expression.
%
%   Example:
%   valid_trials = selecttrial(TE,'Hit==1|FalseAlarm==1');
%
%   See also FILTERTRIALS and PARSEVALIDTRIALS.

%   Edit log: BH 7/5/12

% Position of logical operators
posAND = findstr(selSTR,'&');   % & (AND)
posOR = findstr(selSTR,'|');    % | (OR)
posANDOR = sort([posAND posOR]);  % any operator

% Separate expression into tokens
posCUT(1,:) = [1 posANDOR+1];   % use the operators as delimiters
posCUT(2,:) = [posANDOR-1 length(selSTR)];
NumTokens = size(posCUT,2);
tokens = cell(1,NumTokens);
ctokens = cell(1,NumTokens);
for i = 1:NumTokens
    tokens{i} = selSTR(posCUT(1,i):posCUT(2,i));   % separate elementary expressions
    ctokens{i} = convert_token(tokens{i});       % convert elementary expressions to 'find' statement (operations will be performed on index sets)
end

% Sort operations
if ~isempty(posANDOR)   % if the expression is already elementary, no need to sort

    % Establish presedence of operations (OR takes precedence over AND)
    junk = regexprep(selSTR(posANDOR),'[&]','1');   % & AND 1
    precedence = regexprep(junk,'[|]','2');         % | OR 2
    
    % Sort according to precedence
    NumOperators = length(posANDOR);
    operations = nan(NumOperators,3);
    for iOPS = 1:NumOperators
        operations(iOPS,:) = [str2double(precedence(iOPS)) iOPS iOPS+1];
    end
    operations = sortrows(operations);      % sort according to precedence
    
    % Create set operation according to precedence
    operationSTR = {'intersect','union'};   % operation is performed on index sets
    for iOPS = 1:size(operations,1)
        oper = operations(iOPS,1);
        tok1 = operations(iOPS,2);
        tok2 = operations(iOPS,3);
        
        % Create set operation tokens
        ctokens{tok1} = sprintf('%s(%s , %s)',operationSTR{oper},ctokens{tok1},ctokens{tok2});
        
        % Replace token2 links with token1
        operations(operations(:,2)==tok2,2) = tok1;
        operations(operations(:,3)==tok2,3) = tok1;
    end
end

% The final expression is in the 1st token
expression_to_execute = ctokens{1};

% If CellID is given instead of events structure, load trial events
if ischar(TE)
   TE = load(cellid2fnames(TE,'TrialEvents2')); %#ok<NASGU>
end

% Execute expression
trials = eval(expression_to_execute);

% -------------------------------------------------------------------------
function outtoken = convert_token(intoken)

% This complicated expression somewhat checks the syntax for validity
outtoken = regexprep(intoken, '\s*(\S+)\s*([<>=~]*[=]*)\s*([-+]*\d*[.]*\d*[e]*\d*)\s*','find(TE.$1 $2  $3)');