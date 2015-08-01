function yn = iscellid(cellid)
%ISCELLID   Check if CellBase cell ID exists.
%   YN = ISCELLID(CELLID) returns true if the referred cell ID exists and
%   false if it doesn't. It checks the existance of a .mat file associated
%   with the cell (see VALIDCELLID) and checks whether the cell ID was
%   stored as a new cell ID.
%
%   See also VALIDCELLID.

%   Edit log: BH 5/4/12

% Search among existing cellids
yn = validcellid(cellid,'list');
yn = logical(yn);

% Check for new cellid
global NEWCELLIDS
if ~yn
    if ~isempty(NEWCELLIDS)
        yn = any(strcmp(NEWCELLIDS,cellid));
    end
end