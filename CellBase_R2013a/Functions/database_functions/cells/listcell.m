function  v = listcell(cellid)
%LISTCELL   List all property values associated with a cell.
%   V = LISTCELL(CELLID) reads TheMatrix (see CellBase documentation) and
%   returns all properties for a cell or multiple cells (CELLID).
%
%   See also FINDCELLPOS.

%   Edit log: BH 5/3/12

% Load CellBase
global CELLIDLIST ANALYSES TheMatrix
if isempty(CELLIDLIST)
    load(getpref('cellbase','fname'));
end

% Find cell(s)
if ischar(cellid)
    pos = findcellpos(cellid);
else
    NumCells = length(cellid);
    pos = nan(1,NumCells);
    for i = 1:NumCells
        pos(i) = findcellpos(cellid(i));
    end
end

% Output
if pos == 0
    disp('LISTCELL: cell not found.');
    return
else
    v = TheMatrix(pos,:);  %#ok<NODEF>
end