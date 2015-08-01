function match = findcellstr(cell_string,pattern)
%FINDCELLSTR   STRCMP for CellBase cell IDs.
%   I = FINDCELLSTR(CELLID,LIST) returns the position of the first matching
%   cell ID in LIST (cell array of strings). It returns 0 for no match.
%
%   See also STRCMP.

%   Edit log: BH 6/22/11

% Find
match = 0;
for i = 1:size(cell_string,1)
    if strcmp(cell_string(i),pattern)
        match = i;
        return
    end
end