function  cellids = selectcell(selSTR,varargin)
%SELECTCELL   Find cells satisfying certain criteria.
%   CELLIDS = SELECTCELL(SELSTR) executes the "logical expression" in
%   SELSTR to find cell IDs matching a set of criteria.
% 
%   The "logical expression" is constructed of criteria based on the values
%   of particular properties (type listtag('property') for a list of all
%   properties in CellBase) connected with logical operators & (AND) and | 
%   (OR). For example '"Rate" > 5 | "Selectivity" > 0.2' will return
%   CellIDs whose "Rate" is either greater than 5 or "Selectivity" is
%   greater than 0.2. "Rate" and "Selectivity" must be existing properties
%   in CellBase (they can be added using ADDANALYSIS or INSERTDATA).
%
%   Operations including string properties are also supported, as in the
%   following example:
%   selstr = ['"PV+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
%       'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
%   PV = selectcell(selstr);
%
%   See also GETVALUE and SELECTTRIAL.

%   Edit log: BH 7/6/12, 8/21/12

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
    tokens{i} = selSTR(posCUT(1,i):posCUT(2,i));
    ctokens{i} = convert_token(tokens{i});
end

% Sort operations
if length(ctokens) > 1     % if the expression is already elementary, no need to sort
    
    % Establish presedence of operations (OR takes precedence over AND)
    junk = regexprep(selSTR(posANDOR),'[&]','1');   % & AND 1
    precedence = regexprep(junk,'[|]','2');         % | OR 2
    
    % Sort according to precedence
    NumOperators = length(posANDOR);
    operations = nan(NumOperators,3);
    for iOPS = 1:length(posANDOR)
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
        ctokens{tok1}=sprintf('%s(%s , %s)',operationSTR{oper},ctokens{tok1},ctokens{tok2});
        
        % Replace token2 links with token1
        operations(operations(:,2)==tok2,2) = tok1;
        operations(operations(:,3)==tok2,3) = tok1;
    end
end

% The final expression is in the 1st token
expression_to_execute = ctokens{1};

% Load CellBase
loadcb

% Execute expression
cellidsPOS = eval(expression_to_execute);
cellids = CELLIDLIST(cellidsPOS);

% -------------------------------------------------------------------------
function outtoken = convert_token(intoken)

% This complicated expression somewhat checks the syntax for validity
% outtoken = regexprep(intoken, '\s*"(.+)"\s*([<>=~][=]*)\s*([-+]*\d[.]*\d*[e]*\d*)\s*','find(getvalue(''$1'') $2  $3)');
% outtoken = regexprep(intoken, '\s*"(.+)"\s*([<>=~][=]*)\s*(([-+]*\d[.]*\d*[e]*\d*)|(\''\w*\''))\s*','find(getvalue(''$1'') $2  $3)');
% outtoken = regexprep(intoken,'\s*"(.+)"\s*(.*)','find(getvalue(''$1'')$2)');
outtoken = regexprep(intoken,'\s*(.*)\s*"(.+)"\s*(.*)','find($1getvalue(''$2'')$3)');