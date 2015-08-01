function list = unique_cell(cellarray)
%UNIQUE_CELL   Unique list of cells.
%   L = UNIQUE_CELL(CELLS) returns unique list of animal-session pairs
%   (N-by-2 cell arrays).
%
%   See also UNIQUE, FINDALLCELLS and ADDNEWCELLS.

%   Edit log: AK 4/04, SPR 12/18/09, BH 3/21/11

% Convert to N-by-1 list
lca = length(cellarray);
junk = cell(lca,1);
for i = 1:lca
  junk{i} = strcat(cellarray{i,:});
end

% Output
[u,pos] = unique(junk);
list = cellarray(pos,:);