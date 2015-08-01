function deletecb(cb_name)
%DELETECB   Delete CellBase.
%   DELETECB(NAME) deletes Matlab preferences for the CellBase with the
%   name specified. The files are not deleted. It forces to choose a new
%   default CellBase if the current one is being deleted.
%
%   See also INITCB and CHOOSECB

%   Edit log: BH 5/30/11

% Check input arguments
error(nargchk(0,1,nargin))
if nargin < 1
    cb_name = input('Which CellBase to delete? ','s');
end

% Find CellBase
crcb = findcb(cb_name);
cellbases = getpref('cellbase','cellbases');

% Delete CellBase
if length(cellbases) < 2
    rmpref('cellbase')    % only one CellBase
    return
else
    cellbases(crcb) = [];
    setpref('cellbase','cellbases',cellbases)
end

% Choose a new CellBases if the current was deleted
if isequal(getpref('cellbase','name'),cb_name)
    disp('You have to choose a new default CellBase.')
    choosecb
end