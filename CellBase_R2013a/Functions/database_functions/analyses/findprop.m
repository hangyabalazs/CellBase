function  c = findprop(property)
%FINDPROP   Find a property in CellBase.
%   POS = FINDPROP(PROPNAME) looks for a property in CellBase that
%   matches PROPNAME and returns the column number in TheMatrix of the
%   matching property.
%
%   See also FINDANALYSIS.

%   Edit log: BH 5/4/12

% Input argument check
if ~ischar(property)
    disp('FINDPROP: wrong argument type.')
    return
end

% Load CellBase
load(getpref('cellbase','fname'),'ANALYSES');

% Find position
[pos0 pos] = findanalysis(property); %#ok<ASGLU>

% Output
if pos == 0
    disp('FINDPROP: No matching property found.');
    return
else
    c = ANALYSES(pos(1)).columns(pos(2));
end