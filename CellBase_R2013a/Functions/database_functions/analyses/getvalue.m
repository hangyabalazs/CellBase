function  v = getvalue(property,varargin)
%GETVALUE   Extract all values associated with a property.
%   V = GETVALUE(PROP) returns a vector of values of the specified
%   property (PROP) for all cellids in CellBase .
%
%   V = GETVALUE(PROP,CELLIDS) returns a vector of values of the specified
%   property (PROP) for each cell in CELLIDS.
%
%   See also FINDANALYSIS.

%   Edit log: BH 8/25/2011

% Load CellBase
global CELLIDLIST ANALYSES TheMatrix
if isempty(CELLIDLIST)
    load(getpref('cellbase','fname'));
end

% Check input arguments
if ~ischar(property)
    disp('GETVALUE: Wrong argument type.')
    v = [];
    return
end

% Get the position of the searched property
[pos0 pos] = findanalysis(property); %#ok<ASGLU>
if pos == 0
    disp('GETVALUE: No matching property found.');
    v = [];
    return
end

% Return the list of values
if nargin > 1
    if iscell(varargin{1})  % values will correspond to cellids without reordering
        posCELLS = cellfun(@(d)find(ismember(CELLIDLIST,d)),varargin{1});
    else
        posCELLS = cellfun(@(d)find(ismember(CELLIDLIST,d)),varargin(1));
    end
else
    posCELLS = 1:length(CELLIDLIST);
end
v = TheMatrix(posCELLS,ANALYSES(pos(1)).columns(pos(2)));

% Convert to double if possible
if iscell(v) && all(cellfun(@(s)isequal(numel(s),1),v)&cellfun(@isnumeric,v))
    v = cell2mat(v);
end