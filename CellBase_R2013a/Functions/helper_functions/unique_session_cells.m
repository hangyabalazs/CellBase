function cells = unique_session_cells(cellids)
%UNIQUE_SESSION_CELLS   Select one cell per session.
%   CELLS = UNIQUE_SESSION_CELLS(CELLIDS) performs a 'unique' operation on
%   session IDs corresponding to the cell IDs in the input argument and
%   returns a new list of cell IDs with only one cell per session.
%
%   See also UNIQUE_CELL.

%   Edit log: AK 07/1, SPR 2010/05, BH 2/28/13

% Input argument check
if nargin < 1
    cellids = listtag('cells');
end
if nargin < 2
    unique_tag_length = 1:12;
end

% Seesion IDs
char_cellids = char(cellids);
char_cellids = char_cellids(:,unique_tag_length);

% Find cell IDs with different session IDs
try
    ids = find(sum(abs(diff(char_cellids))'));
    cells = cellids([1 ids+1]);
catch
    [unique_ids,ids] = unique(char_cellids,'rows');
    cells = cellids(ids);
end